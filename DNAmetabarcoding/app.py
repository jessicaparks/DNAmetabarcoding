#!/usr/bin/env python

import click
import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# blast filtering thresholds
EVALUE = 0.001
IDENTITY = 95
COVERAGE = 90

# set taxonomy ranks to include (in order of increasing specificity)
ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

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


def parse_primers(input, fwd, rev):
    # prep forward primer file (add ^ to the start of each sequence)
    # prep reverse primer file (make reverse compliments and add X to the 3' end)
    fwd_records = []
    rev_records = []
    for record in SeqIO.parse(input, 'fasta'):
        fwd_records.append(
            SeqRecord('^' + record.seq, id=record.id, description='')
        )
        rev_records.append(
            SeqRecord(record.seq.reverse_complement() + 'X',
                      id=record.id, description='')
        )
    SeqIO.write(fwd_records, fwd, 'fasta')
    SeqIO.write(rev_records, rev, 'fasta')


def trim_primers(input, output, adapterless, tooshort, forwardprimers, reverseprimers, samplename):
    # run cutadapt to trim 5' end
    subprocess.check_call([
        'cutadapt',
        '-g', f'file:{forwardprimers}',
        '-e', '0.25', 
        '--untrimmed-output', adapterless,
        '-o', f'{TMP}/{samplename}_temp.fastq',
        input
    ])
    # run cutadapt to trim 3' end
    subprocess.check_call([
        'cutadapt',
        '-a', f'file:{reverseprimers}',
        '-e', '0.25',
        '-m', '50',
        '--too-short-output', tooshort,
        '-o', output,
       	f'{TMP}/{samplename}_temp.fastq'
    ])


def run_dada2(input, output):
    subprocess.check_call(['Rscript', 'dada2.R', input, output])
    return pd.read_csv(output, index_col=0).reset_index()


def dada2_taxonomy(input, output, reference):
    subprocess.check_call(['Rscript', 'dada2_taxonomy.R', input, output, reference])


def write_blast_fasta(data, output):
    records = [
        SeqRecord(Seq(x.sequence), id=str(x.index), description='')
        for x in data.itertuples()
    ]
    SeqIO.write(records, output, 'fasta')


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
    # filter blast results
    blast_header = ['qacc', 'sacc', 'qlen', 'slen', 'pident', 'length', 'qcovs', 'staxid']
    blast_data = pd.read_csv(output, header=None, names=blast_header, sep='\t')
    blast_data = blast_data[blast_data['qcovs']>=COVERAGE].reset_index(drop=True)
    blast_data = pd.merge(
        blast_data.groupby('qacc')[['pident']].max().reset_index(),
        blast_data,
        on=['qacc', 'pident'],
        how='inner'
    )
    return blast_data


def run_taxize(input, output):
    subprocess.check_call(['Rscript', 'taxizedb.R', input, output])
    taxa = pd.read_csv(output)
    taxa = (
        taxa[taxa['rank'].isin(ranks)]
        .pivot(index='taxid', columns='rank', values='name')
        .reset_index()
    )
    return taxa


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


def combine_blast_taxize(blast_data, taxa):
    # merge blast results with taxize information
    merged = (
        pd.merge(
            blast_data[['qacc', 'staxid']],
            taxa,
            left_on='staxid',
            right_on='taxid',
            how='left'
        )
        .drop(['staxid', 'taxid'], axis=1)
        .drop_duplicates()
        .reset_index(drop=True)
    )
    # keep taxa from BLAST only down to the last consistent rank
    final_taxa = pd.DataFrame.from_dict(
        merged.groupby('qacc')
        .apply(lambda x: consistent_taxa(x))
        .to_list()
    )
    return final_taxa


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(exists=False), required=True)
@click.option('--primers', type=click.Path(exists=True), required=True)
@click.option('--taxmethod', type=click.Choice(['BLAST', 'DADA2'], case_sensitive=False), required=True)
@click.option('--taxreference', type=click.Choice(dada2_tax_dbs.keys(), case_sensitive=False))
@click.option('--blastdatabase', default='/gpfs_partners/databases/ncbi/blast/nt/nt')
@click.option('--threads', type=int, default=4, show_default=True)
def main(input, output, primers, taxmethod, taxreference, blastdatabase, threads):

    # set input file base name
    base = os.path.basename(input).replace('.gz', '').replace('.fastq', '')
    
    # set file names for primers and primer trimming
    fwdprimers = f'{TMP}/{base}_forwardprimerstemp.fasta'
    revprimers = f'{TMP}/{base}_reverseprimerstemp.fasta'
    trimmed = f'{TMP}/{base}_trimmed.fastq'
    adapterless = f'{TMP}/{base}_untrimmed.fastq'
    tooshort = f'{TMP}/{base}_tooShort.fastq'
    
    # trim primers
    parse_primers(primers, fwdprimers, revprimers)
    trim_primers(input, trimmed, adapterless, tooshort, fwdprimers, revprimers, base)

    # dada2 asv identification
    asv = f'{TMP}/{base}_asv.csv'
    asv_data = run_dada2(trimmed, asv)

    # taxonomic classification
    # via blast
    if taxmethod == 'BLAST':
        # prep input fasta file for blast
        asv_fasta = f'{TMP}/{base}_asv.fasta'
        write_blast_fasta(asv_data, asv_fasta)

        # run blast
        blast_results = f'{TMP}/{base}_blast.tsv'
        blast_data = run_blast(asv_fasta, blast_results, blastdatabase, threads)
        
        # look up taxa
        taxid_list = f'{TMP}/{base}_taxids.csv'
        blast_data[['staxid']].drop_duplicates().to_csv(taxid_list, header=False, index=False)
        taxize_results = f'{TMP}/{base}_taxa.csv'
        taxa = run_taxize(taxid_list, taxize_results)

        # join blast and taxonomy information
        final_taxa = combine_blast_taxize(blast_data, taxa)

        # join blast taxonomy results with ASV data
        output_data = pd.merge(
            asv_data,
            final_taxa.rename(columns={'qacc': 'index'}),
            on='index',
            how='left'
        ).drop('index', axis=1)

    # via dada2
    if taxmethod == 'DADA2':
        taxa = f'{TMP}/{base}_taxa.csv'
        dada2_taxonomy(asv, taxa, dada2_tax_dbs[taxreference])

        # join dada2 taxonomy results with ASV data
        output_data = pd.merge(
            asv_data,
            pd.read_csv(taxa, index_col=0).reset_index().rename(columns={'index': 'sequence'}),
            on='sequence',
            how='left'
        ).drop('index', axis=1)
        output_data.rename(columns={c: c.lower() for c in output_data.columns}, inplace=True)

    # output table with ASVs, abundance, and taxonomy
    output_data.to_csv(output, index=False)
    print(f'Analysis complete. Results found at {output}')

if __name__ == '__main__':
    main()
