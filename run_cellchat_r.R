# run_cellchat_r.R  — Run real CellChat (Jin et al. 2021) on Kang 2018 PBMC
#
# This compares REAL R CellChat output vs NEUROSCOPE's Python CellChat reimplementation.
# The key question: do the capacity-only pairs also have CellChat score = 0 in real R CellChat?
#
# Usage:
#   Rscript run_cellchat_r.R
#
# Requirements:
#   install.packages("BiocManager")
#   BiocManager::install("CellChat")
#   install.packages(c("Seurat", "Matrix", "R.matlab", "dplyr"))

# Add user library to search path
R_LIB <- "C:/Users/arjun/R_libs"
if (dir.exists(R_LIB)) .libPaths(c(R_LIB, .libPaths()))

suppressWarnings(suppressMessages({
  library(Matrix)
  library(dplyr)
}))

cat("=== Real R CellChat on Kang 2018 PBMC ===\n\n")

# ── Install CellChat if not present ──────────────────────────────────────────
if (!requireNamespace("CellChat", quietly=TRUE)) {
  cat("Installing CellChat...\n")
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("CellChat", ask=FALSE)
}
if (!requireNamespace("R.matlab", quietly=TRUE)) {
  install.packages("R.matlab", repos="https://cran.r-project.org")
}

library(CellChat)
library(R.matlab)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading Kang 2018 PBMC data...\n")

mat_data  <- R.matlab::readMat("converted/kang2018/kang_expression_sparse.mat")
expr_mat  <- mat_data$X  # cells x genes, dense float32
meta      <- read.csv("converted/kang2018/kang_cell_metadata.csv")
gene_df   <- read.csv("converted/kang2018/kang_gene_names.csv")
gene_names <- toupper(gene_df$gene_name)

cat(sprintf("  Expression: %d cells x %d genes\n", nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  Metadata: %d rows\n", nrow(meta)))
cat(sprintf("  Donors: %d\n", length(unique(meta$donor_id))))

# Align if needed
if (nrow(expr_mat) != nrow(meta)) {
  n <- min(nrow(expr_mat), nrow(meta))
  expr_mat <- expr_mat[1:n, ]
  meta <- meta[1:n, ]
}

# Transpose to genes x cells (CellChat format)
expr_t <- t(expr_mat)
rownames(expr_t) <- gene_names
colnames(expr_t) <- paste0("cell_", 1:ncol(expr_t))

# Cell identity (cluster labels)
cell_labels <- meta$CellTypeName
names(cell_labels) <- colnames(expr_t)
cell_labels <- as.factor(cell_labels)
cat(sprintf("  Cell types: %s\n", paste(levels(cell_labels), collapse=", ")))

# ── Create CellChat object ────────────────────────────────────────────────────
cat("\nCreating CellChat object...\n")
cellchat <- createCellChat(
  object   = expr_t,
  meta     = data.frame(labels=cell_labels, row.names=colnames(expr_t)),
  group.by = "labels"
)

# Use human CellChat database
CellChatDB <- CellChatDB.human
cellchat@DB  <- CellChatDB

# ── Preprocessing ─────────────────────────────────────────────────────────────
cat("Preprocessing...\n")
cellchat <- subsetData(cellchat)  # filter to DB genes
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ── Compute communication probabilities ──────────────────────────────────────
cat("Computing communication probabilities (this uses trimean, Jin et al. 2021)...\n")
cellchat <- computeCommunProb(cellchat, type="triMean")
cellchat <- filterCommunication(cellchat, min.cells=10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ── Extract significant interactions ─────────────────────────────────────────
cat("Extracting significant interactions...\n")
df <- subsetCommunication(cellchat)  # all significant L-R pairs
cat(sprintf("  CellChat significant interactions: %d\n", nrow(df)))
cat(sprintf("  Unique L-R pairs: %d\n", nrow(unique(df[c("ligand","receptor")]))))

# ── Compare with NEUROSCOPE capacity-only pairs ──────────────────────────────
cat("\nLoading NEUROSCOPE results...\n")
neuro <- read.csv("results/kang2018/channel_rankings.csv")
neuro_sig <- neuro[neuro$FDRCapacity < 0.05, ]
cat(sprintf("  NEUROSCOPE significant (FDR<0.05): %d pairs\n", nrow(neuro_sig)))

cc_lr_sig <- paste(toupper(df$ligand), toupper(df$receptor), sep="_")
neuro_lr   <- paste(toupper(neuro_sig$Ligand), toupper(neuro_sig$Receptor), sep="_")

both      <- intersect(neuro_lr, cc_lr_sig)
cap_only  <- setdiff(neuro_lr, cc_lr_sig)   # NEUROSCOPE finds, CellChat misses
str_only  <- setdiff(cc_lr_sig, neuro_lr)   # CellChat finds, NEUROSCOPE misses

cat(sprintf("\n2x2 CONTINGENCY (Real R CellChat vs NEUROSCOPE):\n"))
cat(sprintf("  Both methods:       %d\n", length(both)))
cat(sprintf("  Capacity-only:      %d  [NOVEL — real CellChat misses]\n", length(cap_only)))
cat(sprintf("  Strength-only:      %d\n", length(str_only)))

if (length(cap_only) > 0) {
  cat("\nCapacity-only pairs (real CellChat misses, NEUROSCOPE detects):\n")
  for (lr in cap_only) {
    parts <- strsplit(lr, "_")[[1]]
    lig <- parts[1]; rec <- paste(parts[-1], collapse="_")
    row <- neuro_sig[toupper(neuro_sig$Ligand)==lig & toupper(neuro_sig$Receptor)==rec, ]
    cap <- if (nrow(row) > 0) row$CondCapacity[1] else NA
    ci_lo <- if ("BootCI_lo" %in% names(row) && nrow(row)>0) row$BootCI_lo[1] else NA
    ci_hi <- if ("BootCI_hi" %in% names(row) && nrow(row)>0) row$BootCI_hi[1] else NA
    cat(sprintf("  %s -> %s  capacity=%.3f  95CI=[%.3f,%.3f]\n",
                lig, rec, ifelse(is.na(cap),0,cap),
                ifelse(is.na(ci_lo),0,ci_lo), ifelse(is.na(ci_hi),0,ci_hi)))
  }
}

# ── Fisher exact test ─────────────────────────────────────────────────────────
n_both     <- length(both)
n_cap_only <- length(cap_only)
n_str_only <- length(str_only)
n_neither  <- 0  # hard to define without full gene list
ft <- fisher.test(matrix(c(n_both, n_cap_only, n_str_only, max(1,n_neither)), 2, 2))
cat(sprintf("\nFisher exact: OR=%.3f  p=%.4f\n", ft$estimate, ft$p.value))

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create("results/kang2018", showWarnings=FALSE, recursive=TRUE)
write.csv(df, "results/kang2018/cellchat_R_interactions.csv", row.names=FALSE)

summary_lines <- c(
  "Real R CellChat vs NEUROSCOPE — Kang 2018 PBMC",
  paste0("CellChat version: ", packageVersion("CellChat")),
  paste0("CellChat significant: ", nrow(df)),
  paste0("NEUROSCOPE significant (FDR<0.05): ", nrow(neuro_sig)),
  paste0("Both methods: ", n_both),
  paste0("Capacity-only (NOVEL): ", n_cap_only),
  paste0("Strength-only: ", n_str_only),
  paste0("Fisher exact OR: ", round(ft$estimate, 3),
         "  p=", round(ft$p.value, 4))
)
writeLines(summary_lines, "results/kang2018/cellchat_R_summary.txt")

cat("\nSaved: results/kang2018/cellchat_R_interactions.csv\n")
cat("Saved: results/kang2018/cellchat_R_summary.txt\n")
cat("\nDone.\n")
