#!/usr/bin/env Rscript
library(phyloseq)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
ASVfile <- args[1]
Taxafile <- args[2]
taxalevel <- args[3]

ASV_mat <- read.csv(ASVfile)
taxa_mat <- read.csv(Taxafile)

row.names(ASV_mat) <- ASV_mat$sequence
ASV_mat <- ASV_mat %>% select (-sequence)

row.names(taxa_mat) <- taxa_mat$sequence
taxa_mat <- taxa_mat %>% select (-sequence) 


ASV_mat <- as.matrix(ASV_mat)
taxa_mat <- as.matrix(taxa_mat)

ASV = otu_table(ASV_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)

carbom <- phyloseq(ASV, TAX)
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
pdf("abundanceplot.pdf")
plot_bar(carbom, fill = taxalevel)
dev.off()
