#!/usr/bin/env Rscript

library(dada2, quietly=TRUE)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]
ref <- args[3]

seqtab <- read.csv(file)

# assign the taxonomy
taxa <- assignTaxonomy(
    seqtab,
    ref,
    multithread = TRUE,
    tryRC = TRUE
)

# write taxonomy to file
dir.create(dirname(output), showWarnings=FALSE, recursive=TRUE)
write.csv(taxa, output)
