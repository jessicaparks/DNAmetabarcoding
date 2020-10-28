#!/usr/bin/env python

import click
import os
import shutil
import subprocess
import pandas as pd
from Bio.Seq import Seq
from pytaxize import ncbi

TMP = '/share/trnL_blast/tmp'
os.makedirs(TMP, exist_ok=True)

# blast filtering thresholds
EVALUE = 0.001
IDENTITY = 95
COVERAGE = 90


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


def run_blast(input, output, database):
    subprocess.check_call([
        'blastn',
        '-db', database,
        '-query', input,
        '-max_target_seqs', '10',
        '-evalue', str(EVALUE),
        '-perc_identity', str(IDENTITY),
        '-num_threads', '8',
        '-outfmt', '6 qacc sacc qlen slen pident length qcovs staxid',
        '-out', output
    ])


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(exists=False), required=True)
@click.option('--primers', type=click.Path(exists=True), required=True)
@click.option('--taxmethod', type=click.Choice(['BLAST', 'DADA2'], case_sensitive=False), required=True)
@click.option('--taxreference', type=click.Choice(['GTDB', 'UNITE_fungi', 'UNITE_eukaryote'], case_sensitive=False))
@click.option('--blastdatabase', default='/gpfs_partners/databases/ncbi/blast/nt/nt')
@click.option('--entrezkey', type=click.Path(exists=True))
def main(input, output, primers, taxmethod, taxreference, entrezkey, blastdatabase):

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
        #read entrez key
        with open(entrezkey) as f:
            os.environ['ENTREZ_KEY'] = f.read().strip()

        #prep input files
        asv_fasta = f'{TMP}/{base}_asv.fasta'
        blast_results = f'{TMP}/{base}_blast.tsv'
        fh = open(asv, "r")
        sequences = open(asv_fasta, "w")
        for line in fh:
            if line.startswith("\"\","):
                pass
            else:
                data = line.replace("\"","")
                data = data.split(",")
                sequences.write(">" + data[0] + "\n" + data[1] + "\n")
        sequences.close()
        fh.close()
        
        #run blast
        run_blast(asv_fasta, blast_results, blastdatabase)

        #filter blast results
        blast_header = ['qacc', 'sacc', 'qlen', 'slen', 'pident', 'length', 'qcovs', 'staxid']
        blast_data = pd.read_csv(blast_results, header=None, names=blast_header, sep='\t')
        blast_data = blast_data[blast_data['qcovs']>=COVERAGE].reset_index(drop=True)
        top_blast_data = pd.merge(
            blast_data.groupby('qacc')[['pident']].max().reset_index(),
            blast_data,
            on=['qacc', 'pident'],
            how='inner'
        )
        top_blast_data.to_csv('test_blast_filtered_results.csv', index=False)
        
        #look up taxa
        ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        staxids = top_blast_data['staxid'].unique().tolist()
        taxa = (
            pd.DataFrame(
                [[x[0],y['Rank'],y['ScientificName']]
                  for x in ncbi.hierarchy(ids=staxids).items()
                  for y in x[1] if y['Rank'] in ranks],
                columns=['staxid', 'rank', 'scientificname']
            )
            .pivot(index='staxid', columns='rank', values='scientificname')
            .reset_index()
        )
        top_blast_data = (
            pd.merge(
                top_blast_data[['qacc', 'staxid']],
                taxa,
                on='staxid',
                how='left'
            )
            .drop('staxid', axis=1)
            .drop_duplicates()
            .reset_index(drop=True)
        )

        #keep taxa from BLAST only down to the last consistent rank
        def consistent_taxa(x):
            new_taxa = {'qacc': x['qacc'].unique().tolist()[0]}
            stop = 0
            for r in ranks:
                r_taxa = x[r].dropna().unique().tolist()
                if len(r_taxa) == 1:
                    new_taxa[r] = r_taxa[0]
                else:
                    stop += 1
                if stop > 0:
                    new_taxa[r] = None
            return(new_taxa)
        final_taxa = pd.DataFrame.from_dict(
            top_blast_data.groupby('qacc')
            .apply(lambda x: consistent_taxa(x))
            .to_list()
        )
        output_data = pd.merge(
            pd.read_csv(asv, index_col=0).reset_index(),
            final_taxa.rename(columns={'qacc': 'index'}),
            on='index',
            how='left'
        )

    # via dada2
    if taxmethod == 'DADA2':
        refs = {
            'GTDB': '/share/trnL_blast/dada2-reference/GTDB.fasta',
            'UNITE_fungi': '/share/trnL_blast/dada2-reference/UNITE_fungi.fasta',
            'UNITE_eukaryote': '/share/trnL_blast/dada2-reference/UNITE_eukaryote.fasta'
        }
        taxa = f'{TMP}/{base}_taxa.csv'
        dada2_taxonomy(asv, taxa, refs[taxreference])
        output_data = pd.merge(
            pd.read_csv(asv, index_col=0).reset_index(),
            pd.read_csv(taxa, index_col=0).reset_index().rename(columns={'index': 'sequence'}),
            on='sequence',
            how='left'
        )
        output_data.rename(columns={c: c.lower() for c in output_data.columns}, inplace=True)

    # output table
    output_data.to_csv(output, index=False)

if __name__ == '__main__':
    main()
