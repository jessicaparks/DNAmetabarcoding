#!/usr/bin/env python

import click
import os
import pandas as pd
import subprocess
from Bio.Seq import Seq

# blast filtering thresholds
EVALUE = 0.001
IDENTITY = 95
COVERAGE = 90

# dada2 taxonomy database paths
base_dada2_tax = '/usr/local/usrapps/trnL_blast/dada2_taxonomy'
dada2_tax_dbs = {
    'GTDB': os.path.join(base_dada2_tax, 'GTDB.fasta'),
    'UNITE_fungi': os.path.join(base_dada2_tax, 'UNITE_fungi.fasta'),
    'UNITE_eukaryote': os.path.join(base_dada2_tax, 'UNITE_eukaryote.fasta')
}

# set temporary	directory for intermediate files
user = os.environ['USER']
TMP = f'/share/trnL_blast/{user}/tmp/dnametabarcoding'
os.makedirs(TMP, exist_ok=True)


def trim_primers(input, output, adapterless, tooShort, forwardprimers, reverseprimers):
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


def run_blast(input, output, database, threads):
    subprocess.check_call([
        'blastn',
        '-db', database,
        '-query', input,
        '-max_target_seqs', '10',
        '-evalue', str(EVALUE),
        '-perc_identity', str(IDENTITY),
        '-num_threads', str(threads),
        '-outfmt', '6 qacc sacc qlen slen pident length qcovs staxid',
        '-out', output
    ])


def run_taxize(input, output):
    subprocess.check_call(['Rscript', 'taxizedb.R', input, output])


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(exists=False), required=True)
@click.option('--primers', type=click.Path(exists=True), required=True)
@click.option('--taxmethod', type=click.Choice(['BLAST', 'DADA2'], case_sensitive=False), required=True)
@click.option('--taxreference', type=click.Choice(dada2_tax_dbs.keys(), case_sensitive=False))
@click.option('--blastdatabase', default='/gpfs_partners/databases/ncbi/blast/nt/nt')
@click.option('--threads', type=int, default=4, show_default=True)
def main(input, output, primers, taxmethod, taxreference, blastdatabase, threads):

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
            if line.startswith("\"\","):
                pass
            else:
                data = line.replace("\"","")
                data = data.split(",")
                sequences.write(">" + data[0] + "\n" + data[1] + "\n")
        sequences.close()
        fh.close()
        
        #run blast
        run_blast(asv_fasta, blast_results, blastdatabase, threads)

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
        
        #look up taxa
        ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        taxid_list = f'{TMP}/{base}_taxids.csv'
        top_blast_data[['staxid']].drop_duplicates().to_csv(taxid_list, header=False, index=False)
        taxize_results = f'{TMP}/{base}_taxa.csv'
        run_taxize(taxid_list, taxize_results)
        taxa = pd.read_csv(taxize_results)
        taxa = (
            taxa[taxa['rank'].isin(ranks)]
            .pivot(index='taxid', columns='rank', values='name')
            .reset_index()
        )
        top_blast_data = (
            pd.merge(
                top_blast_data[['qacc', 'staxid']],
                taxa,
                left_on='staxid',
                right_on='taxid',
                how='left'
            )
            .drop(['staxid', 'taxid'], axis=1)
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
        taxa = f'{TMP}/{base}_taxa.csv'
        dada2_taxonomy(asv, taxa, dada2_tax_dbs[taxreference])
        output_data = pd.merge(
            pd.read_csv(asv, index_col=0).reset_index(),
            pd.read_csv(taxa, index_col=0).reset_index().rename(columns={'index': 'sequence'}),
            on='sequence',
            how='left'
        )
        output_data.rename(columns={c: c.lower() for c in output_data.columns}, inplace=True)

    # output table with ASVs, abundance, and taxonomy
    output_data.to_csv(output, index=False)
    print(f'Results at {output}')

if __name__ == '__main__':
    main()
