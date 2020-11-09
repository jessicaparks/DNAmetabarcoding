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

png("abundanceplot.png")
if(taxalevel == "kingdom"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=kingdom, fill=kingdom), stat="identity", position="stack")
} else if(taxalevel == "phylum"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
} else if(taxalevel == "class"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=class, fill=class), stat="identity", position="stack")
} else if(taxalevel == "order"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=order, fill=order), stat="identity", position="stack")
} else if(taxalevel == "family"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=family, fill=family), stat="identity", position="stack")
} else if(taxalevel == "genus"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")
} else if(taxalevel == "species"){
  plot_bar(carbom, fill = taxalevel) + geom_bar(aes(color=species, fill=species), stat="identity", position="stack")
}
dev.off()
