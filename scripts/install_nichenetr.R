# Install nichenetr from GitHub (saeyslab/nichenetr).
# Installs into your default user library.
# Run: Rscript install_nichenetr.R
user_lib <- Sys.getenv("R_LIBS_USER", unset=file.path(path.expand("~"), "R", "library"))
dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(user_lib, .libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'), timeout=300)
cat('Installing deps...\n')
install.packages(c('RcppArmadillo','remotes','dplyr','ggplot2','tidyr'), lib=user_lib, quiet=TRUE)
library(remotes, lib.loc=user_lib)
remotes::install_github('saeyslab/nichenetr', lib=user_lib, upgrade='never')
cat('DONE\n')