#!/usr/bin/env Rscript

# set HOME environment variable to group directory because
# the taxize database will be downloaded to this location
Sys.setenv(HOME = "/usr/local/usrapps/trnL_blast")

# load taxizedb library
library(taxizedb, quietly=TRUE)

# download the NCBI taxize database
path <- db_download_ncbi(overwrite=TRUE, verbose=TRUE)
cat("NCBI taxize db downloaded to", path, "\n")
