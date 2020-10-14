#!/usr/bin/env Rscript

library(dada2, quietly=TRUE)
suppressPackageStartupMessages(library(DECIPHER, quietly=TRUE))

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]
ref <- args[3]

seqtab <- read.csv(file)

# create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab))

# change to the path of the training set
load(ref)

# assign the taxonomy
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Convert the output object of class "Taxa" to a matrix analogous to the output from DADA2's assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)

# write taxonomy to file
dir.create(dirname(output), showWarnings=FALSE, recursive=TRUE)
write.csv(taxid, output)
