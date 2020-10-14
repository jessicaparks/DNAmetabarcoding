#!/usr/bin/env python

import click
import os
import shutil
import subprocess

TMP = '/share/trnL_blast/tmp'
os.makedirs(TMP, exist_ok=True)


def trim_primers(input, primers, output, adapterless, tooShort):
    shutil.copy(input, output)
    #prep forward primer file(add ^ to the start of each sequence)
    
    #run cutAdapt to trim 5' end
    #cutadapt5Prime = "cutadapt.exe -g file:" + primers + " -e=0.25 --untrimmed-output " + adapterless + " -o temp.fastq " + input
    #os.system(cutadapt5Prime)

    #prep reverse primer file(make reverse compliments and add X to the 3' end)

    #run cutAdapt to trim 3' end
    #cutadapt3Prime = "cutadapt.exe -a file:" + reverseprimers + " -e=0.25 -m=50 --too-short-output " + tooShort + " -o  " + output + "temp.fastq"
    #os.system(cutadapt3Prime)



def run_dada2(input, output):
    subprocess.check_call(['Rscript', 'dada2.R', input, output])


def dada2_taxonomy(input, output, reference):
    subprocess.check_call(['Rscript', 'dada2_taxonomy.R', input, output, reference])


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('--primers', type=click.Path(exists=True), required=True)
@click.option('--taxmethod', type=click.Choice(['DADA2', 'BLAST'], case_sensitive=False), required=True)
@click.option('--taxreference', type=click.Choice(['GTDB', 'UNITE_fungi', 'UNITE_eukaryote'], case_sensitive=False))
def main(input, primers, taxmethod, taxreference):

    # trim primers
    base = os.path.basename(input).replace('.fastq.gz', '')
    trimmed = f'{TMP}/{base}_trimmed.fastq'
    adapterless = f'{TMP}/{base}_untrimmed.fastq'
    tooShort = trimmed = f'{TMP}/{base}_tooShort.fastq'
    trim_primers(input, primers, trimmed, adapterless, tooShort)

    # dada2 asv identification
    asv = f'{TMP}/{base}_asv.csv'
    run_dada2(trimmed, asv)

    # taxonomic classification
    # via dada2
    if taxmethod == 'DADA2':
        refs = {
            'GTDB': '/share/trnL_blast/dada2-reference/GTDB.fasta',
            'UNITE_fungi': '/share/trnL_blast/dada2-reference/UNITE_fungi.fasta',
            'UNITE_eukaryote': '/share/trnL_blast/dada2-reference/UNITE_eukaryote.fasta'
        }
        taxa = f'{TMP}/{base}_taxa.csv'
        dada2_taxonomy(asv, taxa, refs[taxreference])

    # via blast
    if taxmethod == 'BLAST':
        pass

    # filter data

    # output table


if __name__ == '__main__':
    main()
