# compare_nichenet.R  —  NicheNet vs NEUROSCOPE comparison
#
# Runs NicheNet on Kang 2018 and COVID Atlas to benchmark against NEUROSCOPE.
# NicheNet uses prior-knowledge networks (TF-target, PPI, signaling) to rank
# which ligands in sender cells most likely drive gene expression changes in
# receiver cells. Unlike CellChat (trimean) or NEUROSCOPE (KSG MI), NicheNet
# uses gene-set enrichment of downstream targets as the readout.
#
# Install once:
#   install.packages("BiocManager")
#   BiocManager::install("nichenetr")
#   install.packages(c("tidyverse","Seurat"))
#
# Run: Rscript compare_nichenet.R

suppressPackageStartupMessages({
  library(nichenetr)
  library(tidyverse)
})

cat("=== NicheNet vs NEUROSCOPE comparison ===\n\n")

# ── Load NicheNet prior networks ──────────────────────────────────────────────
cat("Loading NicheNet prior networks...\n")
ligand_target_matrix   <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
lr_network             <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
weighted_networks      <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

# ── Helper: run NicheNet for one dataset ──────────────────────────────────────
run_nichenet <- function(expr_csv, meta_csv, neuroscope_csv, dataset_name,
                         sender_celltypes, receiver_celltype,
                         condition_col = "condition",
                         condition_val = NULL) {

  cat(sprintf("\n--- %s ---\n", dataset_name))

  expr <- read.csv(expr_csv, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_csv)

  # Build Seurat-like expression (genes x cells)
  cat(sprintf("  Expression: %d cells x %d genes\n", nrow(expr), ncol(expr)))

  # Filter to condition if specified
  if (!is.null(condition_val) && condition_col %in% colnames(meta)) {
    keep <- meta[[condition_col]] == condition_val
    meta <- meta[keep, ]
    expr <- expr[keep, ]
  }

  cluster_col <- if ("cell_type" %in% colnames(meta)) "cell_type" else "Cluster"
  cat(sprintf("  Cell types: %s\n", paste(unique(meta[[cluster_col]]), collapse=", ")))

  # Receiver cells: find DE genes (genes upregulated in receiver vs other clusters)
  receiver_mask <- meta[[cluster_col]] %in% receiver_celltype
  sender_mask   <- meta[[cluster_col]] %in% sender_celltypes

  if (sum(receiver_mask) < 10 || sum(sender_mask) < 10) {
    cat(sprintf("  Skipping: too few cells in receiver (%d) or senders (%d)\n",
                sum(receiver_mask), sum(sender_mask)))
    return(NULL)
  }

  # DE genes in receiver vs other clusters
  expr_mat <- t(as.matrix(expr))   # genes x cells
  receiver_expr <- expr_mat[, receiver_mask, drop=FALSE]
  other_expr    <- expr_mat[, !receiver_mask, drop=FALSE]

  # Wilcoxon to find geneset of interest (receiver upregulated genes)
  de_pval <- apply(expr_mat, 1, function(row) {
    tryCatch(wilcox.test(row[receiver_mask], row[!receiver_mask],
                         alternative="greater")$p.value,
             error=function(e) 1)
  })
  de_genes <- names(sort(de_pval)[1:min(500, length(de_pval))])
  cat(sprintf("  DE genes in receiver: %d\n", length(de_genes)))

  # NicheNet ligand activity
  ligands_present <- intersect(rownames(ligand_target_matrix),
                                colnames(expr_mat)[sender_mask])
  cat(sprintf("  Ligands in sender cells: %d\n", length(ligands_present)))

  if (length(ligands_present) == 0 || length(de_genes) == 0) {
    cat("  Skipping: no overlap with NicheNet network.\n")
    return(NULL)
  }

  ligand_activities <- predict_ligand_activities(
    geneset            = de_genes,
    background_expressed_genes = rownames(ligand_target_matrix)[
      rowMeans(expr_mat > 0) > 0.05],
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands  = ligands_present
  )

  ligand_activities <- ligand_activities %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(rank_nichenet = row_number())

  # Load NEUROSCOPE results
  neuroscope_res <- read.csv(neuroscope_csv)
  neuro_ligs <- toupper(neuroscope_res$Ligand)

  # Find which NicheNet top ligands match NEUROSCOPE
  top_n   <- min(20, nrow(ligand_activities))
  top_nn  <- toupper(ligand_activities$test_ligand[1:top_n])

  neuro_only <- setdiff(neuro_ligs, top_nn)
  nn_only    <- setdiff(top_nn, neuro_ligs)
  both       <- intersect(neuro_ligs, top_nn)

  cat(sprintf("\n  NicheNet top %d ligands vs NEUROSCOPE:\n", top_n))
  cat(sprintf("    Both:           %d\n", length(both)))
  cat(sprintf("    NEUROSCOPE-only:%d\n", length(neuro_only)))
  cat(sprintf("    NicheNet-only:  %d\n", length(nn_only)))

  # Build comparison table
  comp <- data.frame(
    Ligand          = neuroscope_res$Ligand,
    Receptor        = if("Receptor" %in% colnames(neuroscope_res)) neuroscope_res$Receptor else NA,
    NEUROSCOPE_Cap  = if("CondCapacity" %in% colnames(neuroscope_res)) neuroscope_res$CondCapacity else NA,
    NEUROSCOPE_FDR  = if("FDRCapacity" %in% colnames(neuroscope_res)) neuroscope_res$FDRCapacity else NA,
    In_NicheNet_Top20 = toupper(neuroscope_res$Ligand) %in% top_nn,
    NicheNet_AUPR   = ligand_activities$aupr_corrected[
      match(toupper(neuroscope_res$Ligand), toupper(ligand_activities$test_ligand))],
    NicheNet_Rank   = ligand_activities$rank_nichenet[
      match(toupper(neuroscope_res$Ligand), toupper(ligand_activities$test_ligand))]
  )

  return(list(
    ligand_activities = ligand_activities,
    comparison        = comp,
    n_both            = length(both),
    n_neuro_only      = length(neuro_only),
    n_nn_only         = length(nn_only)
  ))
}

# ── Kang 2018 ─────────────────────────────────────────────────────────────────
cat("\nRunning NicheNet on Kang 2018...\n")

# Need expression as genes x cells CSV — generate from the NEUROSCOPE channel_rankings
# as a proxy (NicheNet needs raw expression; we use the channel_rankings for ligand list)
kang_out <- tryCatch({
  run_nichenet(
    expr_csv        = "converted/kang2018/kang_expression_dense.csv",  # if available
    meta_csv        = "converted/kang2018/kang_cell_metadata.csv",
    neuroscope_csv  = "results/kang2018/channel_rankings.csv",
    dataset_name    = "Kang 2018 IFN-beta PBMC",
    sender_celltypes = c("CD14+ Monocyte", "DC", "CD4 T"),
    receiver_celltype = c("CD4 T", "CD8 T", "NK")
  )
}, error = function(e) {
  cat(sprintf("  Kang 2018 NicheNet error: %s\n", e$message))
  NULL
})

if (!is.null(kang_out)) {
  write.csv(kang_out$comparison,
            "results/kang2018/comparison_nichenet_vs_neuroscope.csv",
            row.names = FALSE)
  cat(sprintf("  Saved: results/kang2018/comparison_nichenet_vs_neuroscope.csv\n"))
}

# ── COVID Atlas ───────────────────────────────────────────────────────────────
cat("\nRunning NicheNet on COVID Atlas...\n")

covid_out <- tryCatch({
  run_nichenet(
    expr_csv        = "converted/covid_atlas/covid_expression_dense.csv",
    meta_csv        = "converted/covid_atlas/covid_cell_metadata.csv",
    neuroscope_csv  = "results/covid_atlas/channel_rankings.csv",
    dataset_name    = "COVID Atlas (Stephenson 2021)",
    sender_celltypes = c("CD14+ Monocyte", "CD16+ Monocyte", "DC"),
    receiver_celltype = c("CD4 T", "NK", "B")
  )
}, error = function(e) {
  cat(sprintf("  COVID Atlas NicheNet error: %s\n", e$message))
  NULL
})

if (!is.null(covid_out)) {
  write.csv(covid_out$comparison,
            "results/covid_atlas/comparison_nichenet_vs_neuroscope.csv",
            row.names = FALSE)
  cat(sprintf("  Saved: results/covid_atlas/comparison_nichenet_vs_neuroscope.csv\n"))
}

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== NicheNet comparison summary ===\n")
cat("NicheNet uses prior-knowledge networks to rank ligands by downstream\n")
cat("gene-set enrichment in receiver cells — orthogonal to KSG MI.\n\n")
cat("Key distinction:\n")
cat("  NicheNet: ligand -> TF -> target gene enrichment (knowledge-driven)\n")
cat("  CellChat: geometric mean co-expression (data-driven, zero-floor)\n")
cat("  NEUROSCOPE: KSG mutual information channel capacity (data-driven, no zero-floor)\n\n")
cat("Overlap with NEUROSCOPE is expected to be partial:\n")
cat("  - Pairs found by both: high-confidence (two independent methods agree)\n")
cat("  - NEUROSCOPE-only: information-theoretically active, not in prior networks\n")
cat("  - NicheNet-only: prior-knowledge supported, but may lack statistical power\n")
cat("    in sparse single-cell data\n")
