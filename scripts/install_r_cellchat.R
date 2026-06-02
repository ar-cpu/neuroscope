# install_r_cellchat.R — Install all R dependencies for CellChat comparison
#
# Run this once:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" install_r_cellchat.R

R_LIB <- Sys.getenv("R_LIBS_USER", unset=file.path(path.expand("~"), "R", "library"))
dir.create(R_LIB, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(R_LIB, .libPaths()))
cat("Library:", R_LIB, "\n\n")

options(repos="https://cloud.r-project.org")

# CRAN packages
cran_pkgs <- c("BiocManager", "Matrix", "R.matlab", "dplyr", "ggplot2",
               "devtools", "remotes", "NMF", "igraph", "sna",
               "stringr", "reshape2", "patchwork", "cowplot")

cat("Installing CRAN packages...\n")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    install.packages(pkg, lib=R_LIB, quiet=TRUE)
    cat("  Installed:", pkg, "\n")
  } else {
    cat("  Already present:", pkg, "\n")
  }
}

# Bioconductor packages (CellChat dependencies)
library(BiocManager, lib.loc=R_LIB)
bioc_pkgs <- c("BiocNeighbors", "SingleCellExperiment", "MAST")
cat("\nInstalling Bioconductor packages...\n")
BiocManager::install(bioc_pkgs, lib=R_LIB, ask=FALSE, update=FALSE)

# CellChat from GitHub
library(remotes, lib.loc=R_LIB)
cat("\nInstalling CellChat from GitHub (sqjin/CellChat)...\n")
remotes::install_github("sqjin/CellChat", lib=R_LIB, upgrade="never")

cat("\nVerifying CellChat install...\n")
if (requireNamespace("CellChat", quietly=TRUE, lib.loc=R_LIB)) {
  cat("CellChat", as.character(packageVersion("CellChat", lib.loc=R_LIB)), "installed successfully\n")
} else {
  cat("ERROR: CellChat not found after install\n")
}

cat("\nAll dependencies installed. Run: Rscript run_cellchat_r.R\n")
