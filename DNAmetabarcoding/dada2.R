#!/usr/bin/env Rscript

library(dada2)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]
samplename <- gsub(pattern='.fastq.gz','',basename(file))

# filter and trim the reads (DADA2 requires no Ns)
filterAndTrim(
    fwd = file,
    filt = 'filtered.fastq.gz',
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    minLen = 50,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
)

# learn the error rates
err <- learnErrors(
    fls = 'filtered.fastq.gz',
    multithread = TRUE
)

# dereplicate identical reads
derep <- derepFastq(fls = 'filtered.fastq.gz')

# apply sample inference
dadaresult <- dada(
    derep = derep,
    err = err,
    multithread = TRUE
)

# create sequence table and remove chimeras
seqtab <- makeSequenceTable(
    samples = list(samplename=dadaresult)
)
seqtab.nochim <- removeBimeraDenovo(
    unqs = seqtab,
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
)

# write sequence table to file
dir.create(dirname(output), showWarnings=FALSE, recursive=TRUE)
write.csv(seqtab.nochim, output)
