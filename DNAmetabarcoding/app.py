#!/usr/bin/env python

import click
import os
import shutil
import subprocess
from Bio.Seq import Seq

TMP = '/share/trnL_blast/tmp'
os.makedirs(TMP, exist_ok=True)


def trim_primers(input, output, adapterless, tooShort, forwardprimers, reverseprimers):
    #shutil.copy(input, output)
    
    #run cutAdapt to trim 5' end
    cutadapt5Prime = "cutadapt -g file:" + forwardprimers + " -e=0.25 --untrimmed-output " + adapterless + f" -o {TMP}/temp.fastq " + input
    os.system(cutadapt5Prime)

    #run cutAdapt to trim 3' end
    cutadapt3Prime = "cutadapt -a file:" + reverseprimers + " -e=0.25 -m=50 --too-short-output " + tooShort + " -o  " + output + f" {TMP}/temp.fastq"
    os.system(cutadapt3Prime)

def run_dada2(input, output):
    subprocess.check_call(['Rscript', 'dada2.R', input, output])


def dada2_taxonomy(input, output, reference):
    subprocess.check_call(['Rscript', 'dada2_taxonomy.R', input, output, reference])


def idtaxa_taxonomy(input, output, reference):
    subprocess.check_call(['Rscript', 'idtaxa_taxonomy.R', input, output, reference])


def run_blast(input, output):
    subprocess.check_call([
        'blastn',
        '-db', '/gpfs_partners/databases/ncbi/blast/nt/nt',
        '-query', input,
        '-max_target_seqs', '50',
        '-num_threads', '8',
        '-outfmt', '6 qacc sacc qlen slen pident length qcovs staxid ssciname',
        '-out', output
    ])


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('--primers', type=click.Path(exists=True), required=True)
@click.option('--taxmethod', type=click.Choice(['BLAST', 'DADA2', 'IDTAXA'], case_sensitive=False), required=True)
@click.option('--taxreference', type=click.Choice(['GTDB', 'UNITE', 'UNITE_fungi', 'UNITE_eukaryote'], case_sensitive=False))
def main(input, primers, taxmethod, taxreference):

    # set file path name
    base = os.path.basename(input).replace('.fastq.gz', '')
    fh = open(primers, "r")
    
    #prep forward primer file(add ^ to the start of each sequence)
    fhforward = open(f"{TMP}/forwardprimerstemp.fasta", "w")
    for line in fh:
        if line.startswith(">") or line == "\n":
            fhforward.write(line)
        else:
            fhforward.write("^" + line)
    fhforward.close()
    fh.close()

    #prep reverse primer file(make reverse compliments and add X to the 3' end)
    fh = open(primers, "r")
    fhreverse = open(f"{TMP}/reverseprimerstemp.fasta", "w")
    for line in fh:
        if line.startswith(">") or line == "\n":
            fhreverse.write(line)
        else:
            seq = Seq(line.strip("\n"))
            reverseSeq = seq.reverse_complement()
            fhreverse.write(str(reverseSeq + "X" + "\n"))
    fhreverse.close()
    fh.close()
    
    #set file names for trimming
    trimmed = f'{TMP}/{base}_trimmed.fastq'
    adapterless = f'{TMP}/{base}_untrimmed.fastq'
    tooShort = f'{TMP}/{base}_tooShort.fastq'
    
    #trim primers
    trim_primers(input, trimmed, adapterless, tooShort, f"{TMP}/forwardprimerstemp.fasta", f"{TMP}/reverseprimerstemp.fasta")

    # dada2 asv identification
    asv = f'{TMP}/{base}_asv.csv'
    run_dada2(trimmed, asv)

    # taxonomic classification
    # via blast
    if taxmethod == 'BLAST':
        #prep input files
        asv_fasta = f'{TMP}/{base}_asv.fasta'
        blast_results = f'{TMP}/{base}_blast.tsv'
        fh = open(asv, "r")
        sequences = open(asv_fasta, "w")
        for line in fh:
            if line.startswith("\"\",):
                pass
            else:
                data = line.strip("\"")
                data = data.split(",")
                sequences.write(">" + data[0] + "\n" + data[1] + "\n")
        sequences.close()
        fh.close()
        
        #run blast
        run_blast(asv_fasta, blast_results)

    # via dada2
    if taxmethod == 'DADA2':
        refs = {
            'GTDB': '/share/trnL_blast/dada2-reference/GTDB.fasta',
            'UNITE_fungi': '/share/trnL_blast/dada2-reference/UNITE_fungi.fasta',
            'UNITE_eukaryote': '/share/trnL_blast/dada2-reference/UNITE_eukaryote.fasta'
        }
        taxa = f'{TMP}/{base}_taxa.csv'
        dada2_taxonomy(asv, taxa, refs[taxreference])

    # via idtaxa
    if taxmethod == 'IDTAXA':
        refs = {
            'GTDB': '/share/trnL_blast/idtaxa-reference/GTDB_r95-mod_August2020.RData',
            'UNITE': '/share/trnL_blast/idtaxa-reference/UNITE_v2020_February2020.RData'
        }
        taxa = f'{TMP}/{base}_taxa.csv'
        idtaxa_taxonomy(asv, taxa, refs[taxreference])

    # filter data

    # output table


if __name__ == '__main__':
    main()
