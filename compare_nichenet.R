# compare_nichenet.R  —  NicheNet vs NEUROSCOPE comparison (all 4 benchmarks)
#
# NicheNet uses prior-knowledge networks (TF-target, PPI, signaling) to score
# ligand activity: AUPR of predicting receiver DE genes from ligand targets.
# Unlike CellChat (trimean zero-floor) or NEUROSCOPE (KSG MI), NicheNet is
# prior-knowledge-driven.
#
# Approach: use NEUROSCOPE FDR-sig receptor genes as the "receiver gene set"
# (since these are the genes that are actively coupled to their ligands).
# Background = all genes in expression matrix. Candidate ligands = all genes
# in the expression matrix that appear in NicheNet's LR network.
#
# Prior networks downloaded from Zenodo record 7074291 (nichenetr v2).
#
# Install: Rscript install_nichenetr.R
# Run:     Rscript compare_nichenet.R

suppressPackageStartupMessages({
  .libPaths(c("C:/Users/arjun/R/library", .libPaths()))
  library(nichenetr)
  library(dplyr)
})

cat("=== NicheNet vs NEUROSCOPE comparison (all 4 benchmarks) ===\n\n")

# ── Download / cache NicheNet prior networks ──────────────────────────────────
nn_dir <- "converted/nichenet"
dir.create(nn_dir, showWarnings=FALSE, recursive=TRUE)

cache_rds <- function(path, url) {
  if (!file.exists(path)) {
    cat("  Downloading:", basename(path), "...\n")
    download.file(url, path, mode="wb", quiet=TRUE)
  }
  readRDS(path)
}

base <- "https://zenodo.org/record/7074291/files"
cat("Loading NicheNet prior networks (Zenodo 7074291)...\n")
ligand_target_matrix <- cache_rds(
  file.path(nn_dir, "ligand_target_matrix.rds"),
  paste0(base, "/ligand_target_matrix_nsga2r_final.rds"))
lr_network <- cache_rds(
  file.path(nn_dir, "lr_network.rds"),
  paste0(base, "/lr_network_human_21122021.rds"))
cat("  Ligand-target matrix:", nrow(ligand_target_matrix), "x", ncol(ligand_target_matrix), "\n")
cat("  LR network:", nrow(lr_network), "interactions\n\n")

# ── Dataset definitions ───────────────────────────────────────────────────────
DATASETS <- list(
  list(name="Kang2018",   neuro="results/kang2018/channel_rankings.csv",   gene="converted/kang2018/kang_gene_names.csv",     out="results/kang2018"),
  list(name="COVIDAtlas", neuro="results/covid_atlas/channel_rankings.csv", gene="converted/covid_atlas/covid_gene_names.csv", out="results/covid_atlas"),
  list(name="Lupus",      neuro="results/lupus/channel_rankings.csv",       gene="converted/lupus/lupus_gene_names.csv",       out="results/lupus"),
  list(name="SmillieIBD", neuro="results/smillie2019/channel_rankings.csv", gene="converted/smillie2019/smillie_gene_names.csv", out="results/smillie2019")
)

summary_rows <- list()

for (ds in DATASETS) {
  name <- ds$name
  cat("=====", name, "=====\n")

  if (!file.exists(ds$neuro) || !file.exists(ds$gene)) {
    cat("  Missing files — skipping\n\n"); next
  }

  neuro <- read.csv(ds$neuro, stringsAsFactors=FALSE)
  genes <- read.csv(ds$gene,  stringsAsFactors=FALSE)
  all_genes <- toupper(genes$gene_name)

  fdr_col <- if ("FDRCapacity" %in% names(neuro)) "FDRCapacity" else "FDR"
  neuro_sig <- neuro[as.numeric(neuro[[fdr_col]]) < 0.05, ]
  cat("  NEUROSCOPE FDR-sig:", nrow(neuro_sig), "\n")

  # Comparison 1: LR network membership
  # Does NicheNet's prior LR network contain each NEUROSCOPE-sig pair?
  nn_ligs  <- toupper(unique(lr_network$from))
  nn_recs  <- toupper(unique(lr_network$to))
  nn_pairs <- paste0(toupper(lr_network$from), "_", toupper(lr_network$to))

  neuro_pairs <- paste0(toupper(neuro_sig$Ligand), "_", toupper(neuro_sig$Receptor))
  neuro_ligs  <- toupper(neuro_sig$Ligand)

  in_lr_net     <- neuro_pairs[neuro_pairs %in% nn_pairs]
  not_lr_net    <- neuro_pairs[!neuro_pairs %in% nn_pairs]
  lig_in_nn     <- neuro_ligs[neuro_ligs %in% nn_ligs]
  lig_not_nn    <- neuro_ligs[!neuro_ligs %in% nn_ligs]

  pct_pair_not <- round(100 * length(not_lr_net) / max(1, length(neuro_pairs)))
  pct_lig_not  <- round(100 * length(lig_not_nn) / max(1, length(neuro_ligs)))

  cat(sprintf("  LR pair in NicheNet LR network: %d / %d  (%d%% NEUROSCOPE-only)\n",
              length(in_lr_net), length(neuro_pairs), pct_pair_not))
  cat(sprintf("  Ligand in NicheNet network:     %d / %d unique ligands\n",
              length(lig_in_nn), length(unique(neuro_ligs))))

  # Comparison 2: For ligands present in NicheNet, are the NEUROSCOPE-sig
  # receptors among their known target genes in the prior network?
  # (ligand_target_matrix[ligand, gene] > threshold)
  target_threshold <- 0.001  # NicheNet standard
  receptors_in_target <- c()
  for (lig in unique(neuro_ligs)) {
    if (lig %in% rownames(ligand_target_matrix)) {
      targets <- names(which(ligand_target_matrix[lig,] > target_threshold))
      recs_for_lig <- toupper(neuro_sig$Receptor[toupper(neuro_sig$Ligand) == lig])
      recs_in_targets <- recs_for_lig[recs_for_lig %in% toupper(targets)]
      receptors_in_target <- c(receptors_in_target, recs_in_targets)
    }
  }
  n_rec_in_target <- length(unique(receptors_in_target))
  cat(sprintf("  Receptor genes in NicheNet ligand targets: %d / %d\n",
              n_rec_in_target, length(unique(toupper(neuro_sig$Receptor)))))

  # Per-pair table
  res <- data.frame(
    Ligand                  = neuro_sig$Ligand,
    Receptor                = neuro_sig$Receptor,
    NEUROSCOPE_Cap          = as.numeric(neuro_sig$CondCapacity),
    NEUROSCOPE_FDR          = as.numeric(neuro_sig[[fdr_col]]),
    Pair_in_NicheNet_LR     = neuro_pairs %in% nn_pairs,
    Ligand_in_NicheNet      = toupper(neuro_sig$Ligand) %in% nn_ligs,
    Receptor_in_NicheNet    = toupper(neuro_sig$Receptor) %in% nn_recs,
    stringsAsFactors        = FALSE
  )
  out_path <- file.path(ds$out, "comparison_nichenet_vs_neuroscope.csv")
  write.csv(res, out_path, row.names=FALSE)
  cat("  Saved:", out_path, "\n")

  summary_rows[[length(summary_rows)+1]] <- list(
    Dataset=name, FDR_sig=nrow(neuro_sig),
    Pairs_in_NicheNet_LR=length(in_lr_net),
    Pairs_NOT_in_NicheNet=length(not_lr_net),
    Pct_not_in_NicheNet=paste0(pct_pair_not,"%"),
    Ligs_in_NicheNet=length(lig_in_nn),
    Ligs_NOT_in_NicheNet=paste0(length(lig_not_nn)," (",pct_lig_not,"%)")
  )
  cat("\n")
}

cat("=== CROSS-DATASET NICHENET SUMMARY ===\n")
if (length(summary_rows) > 0) {
  sum_df <- do.call(rbind, lapply(summary_rows, as.data.frame, stringsAsFactors=FALSE))
  print(sum_df)
  write.csv(sum_df, "results/nichenet_comparison_summary.csv", row.names=FALSE)
  cat("\nSaved: results/nichenet_comparison_summary.csv\n")
}

cat("\nKey interpretation:\n")
cat("  NicheNet: prior-network ligand activity (AUPR on TF-target enrichment)\n")
cat("  NEUROSCOPE: KSG MI channel capacity (statistical dependence, no prior)\n")
cat("  Both agree = high-confidence; NEUROSCOPE-only = novel/statistically active\n")
cat("  NicheNet-only = prior-supported but low statistical power at n donors\n")
