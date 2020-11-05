#!/usr/bin/env Rscript

# set HOME environment variable to group directory because
# the taxize database will be downloaded to this location
Sys.setenv(HOME = "/usr/local/usrapps/trnL_blast")

# load taxizedb library
library(taxizedb, quietly=TRUE)

# download the NCBI taxize database
path <- db_download_ncbi(verbose=TRUE)
cat("NCBI taxize db downloaded to", path, "\n")

# move database
newpath <- "/usr/local/usrapps/trnL_blast/taxizedb/"
c <- file.copy(path, newpath, overwrite=TRUE)
r <- file.remove(path)
cat("NCBI taxize db moved to", paste(newpath, basename(path), sep=""), "\n")
