#!/usr/bin/env python

import click
import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# set blast filtering thresholds
EVALUE = 0.001
IDENTITY = 95
COVERAGE = 90

# set taxonomy ranks to include (in order of increasing specificity)
ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# set dada2 taxonomy database paths
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
    """Prepare the primer sequence for use in trimming.
    A ^ symbol is added to the start of the sequence for use as a forward
    primer, requiring that it is at the start of the read. An X is added
    to the end of the reverse complement sequence for use as a reverse primer,
    requiring that it is at the end of the read and that it can be a partial
    primer sequence.

    Arguments:
    input -- a fasta file containing the primer sequences
    fwd -- output file path for forward primer sequences
    rev -- output file path for reverse complement primer sequences
    """
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


def trim_primers(input, output, primerless, tooshort, forwardprimers,
                 reverseprimers, samplename, cutoff):
    """Trim primer sequences from the sequencing reads.
    The reads without a primer at the start of the read are excluded and
    output to a separate file. The reads that are too short (less than 50
    bases) after trimming are also excluded and output to a separate file.

    Arguments:
    input -- the fastq file containing the sequencing reads
    output -- the file path for the trimmed reads
    primerless -- the file path for the reads that do not have a primer
    tooshort -- the file path for the reads that are too short after trimming
    forwardprimers -- the file path for the forward primers, to be removed
      from the 5' end of the reads
    reverseprimers -- the file path for the reverse complement primers, to be
      removed from the 3' end of the reads
    samplename -- the name for the sample, to be used in naming files
    """
    # run cutadapt to trim 5' end
    subprocess.check_call([
        'cutadapt',
        '-g', f'file:{forwardprimers}',
        '-e', '0.25', 
        '--untrimmed-output', primerless,
        '-o', f'{TMP}/{samplename}_temp.fastq',
        input
    ])
    # run cutadapt to trim 3' end
    subprocess.check_call([
        'cutadapt',
        '-a', f'file:{reverseprimers}',
        '-e', '0.25',
        '-m', str(cutoff),
        '--too-short-output', tooshort,
        '-o', output,
       	f'{TMP}/{samplename}_temp.fastq'
    ])


def run_dada2(input, output, tempdir=TMP):
    """Run DADA2 on the input sequence reads to identify ASVs.

    Arguments:
    input -- the fastq file of sequencing reads
    output -- the file path for the DADA2 results
    tempdir -- directory for temporary, intermediate files

    Returns:
    a dataframe of ASVs with the columns: index, sequence, and abundance
    """
    subprocess.check_call(['Rscript', 'dada2.R', input, output, tempdir])
    return pd.read_csv(output, index_col=0).reset_index()


def dada2_taxonomy(input, output, reference):
    """
    Run DADA2's assignTaxonomy method for taxonomic classification.

    Arguments:
    input -- the file of ASVs output from DADA2
    output -- the file path for the DADA2 taxonomy results
    reference -- the file path for the reference database to use for
      taxonomic assignment

    Returns:
    a dataframe containing the sequence and the taxonomy assignment
    """
    subprocess.check_call([
        'Rscript', 'dada2_taxonomy.R', input, output, reference
    ])
    # read taxonomy results
    taxa = (
        pd.read_csv(output, index_col=0)
        .reset_index()
        .rename(columns={'index': 'sequence'})
    )
    # change column names to lowercase
    taxa.rename(columns={c: c.lower() for c in taxa.columns}, inplace=True)
    # remove the prefix from the taxon names (ie. 'g__' for genus names)
    for c in ranks:
        taxa[c] = taxa[c].str.split('__').str[1]
    return taxa


def write_blast_fasta(data, output):
    """Create a fasta file from a dataframe for use in BLAST.

    Arguments:
    data -- a dataframe containing sequence and index (identifier) columns
      for the ASVs
    output -- the file path for the fasta file
    """
    records = [
        SeqRecord(Seq(x.sequence), id=str(x.index), description='')
        for x in data.itertuples()
    ]
    SeqIO.write(records, output, 'fasta')


def run_blast(input, output, database, threads):
    """Run BLAST and filter the results.
    The BLAST results are filtered by evalue, percent identity, and coverage.
    The result(s) with maximum percent identity are returned for each input
    sequence.

    Arguments:
    input -- the fasta file to be queried with BLAST
    output -- the file path for the BLAST results
    database -- the BLAST database to be searched
    threads -- the threads (or cpus) to use for BLAST

    Returns:
    a dataframe containing the filtered BLAST results
    """
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
    """Run taxizedb to look up the taxa for the input taxids.
    taxizedb uses a local database to look up the NCBI taxonomy.

    Arguments:
    input -- the file containing the list of the taxids
    output -- the file path for the taxizedb results

    Returns:
    a dataframe containing the taxids and corresponding taxa at all ranks
    """
    subprocess.check_call(['Rscript', 'taxizedb.R', input, output])
    taxa = pd.read_csv(output)
    taxa = (
        taxa[taxa['rank'].isin(ranks)]
        .pivot(index='taxid', columns='rank', values='name')
        .reset_index()
    )
    # remove the genus portion that is included in the species name
    taxa['species'] = taxa.apply(
        lambda x:
        x.species.replace(str(x.genus), '').strip(),
        axis=1
    )
    return taxa


def consistent_taxa(x):
    """Determine a consistent taxon at each rank.
    If all results are the same at a given rank, that name is returned for the
    taxon; otherwise, a null value is returned.

    Arguments:
    x -- a dataframe containing the taxa identified for one query sequence

    Returns:
    a dictionary of the consistent taxa names
    """
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
    return new_taxa


def combine_blast_taxize(blast_data, taxa):
    """Combine the BLAST and taxizedb data.
    For each query sequence, a single consistent set of taxon names are
    returned, or a null value for ranks where the taxon names are not
    consistent.

    Arguments:
    blast_data -- the dataframe containing the BLAST results
    taxa -- the dataframe containing the taxa assignments from taxizedb

    Returns:
    a dataframe containing the resulting taxa assignments
    """
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
@click.option('--taxmethod',
              type=click.Choice(['BLAST', 'DADA2'], case_sensitive=False),
              required=True)
@click.option('--taxreference',
              type=click.Choice(dada2_tax_dbs.keys(), case_sensitive=False))
@click.option('--blastdatabase',
              default='/gpfs_partners/databases/ncbi/blast/nt/nt')
@click.option('--threads', type=int, default=4, show_default=True)
def main(input, output, primers, taxmethod, taxreference, blastdatabase, threads, cutoff=50):
    """Identify ASVs and assign taxonomy.
    This is the main function for the app and it's arguments are specified
    from the command line. The output is a CSV file containing the ASVs,
    their abundance, and their taxonomy assignment.

    Arguments:
    input -- fastq file path for sequence reads
    output -- output CSV file path
    primers -- primers fasta file path
    taxmethod -- the method for taxonomy assignment
    taxreference -- the taxonomy reference database to be used if assigning
      taxonomy using DADA2
    blastdatabase -- the path for the BLAST database to be used if assigning
      taxonomy using BLAST
    threads -- the threads/cpus to be used for running multithreaded tasks
    cutoff -- cutoff length for sequences after primer trimming
    """

    # set input file base name and get output directory
    base = os.path.basename(input).replace('.gz', '').replace('.fastq', '')
    outdir = os.path.dirname(output)
    
    # set file names for primers and primer trimming
    fwdprimers = f'{TMP}/{base}_forwardprimerstemp.fasta'
    revprimers = f'{TMP}/{base}_reverseprimerstemp.fasta'
    trimmed = f'{TMP}/{base}_trimmed.fastq'
    primerless = f'{outdir}/{base}_untrimmed.fastq'
    tooshort = f'{outdir}/{base}_tooshort.fastq'
    
    # trim primers
    print('Trimming primers ...')
    parse_primers(primers, fwdprimers, revprimers)
    trim_primers(input, trimmed, primerless, tooshort, fwdprimers, revprimers, base, cutoff)

    # dada2 asv identification
    print('Identifying ASVs with DADA2 ...')
    asv = f'{TMP}/{base}_asv.csv'
    asv_data = run_dada2(trimmed, asv)

    # taxonomic classification
    # via blast
    if taxmethod == 'BLAST':
        # prep input fasta file for blast
        asv_fasta = f'{TMP}/{base}_asv.fasta'
        write_blast_fasta(asv_data, asv_fasta)

        # run blast
        print('Running BLAST ...')
        blast_results = f'{TMP}/{base}_blast.tsv'
        blast_data = run_blast(asv_fasta, blast_results, blastdatabase, threads)
        
        # look up taxa
        print('Running Taxize ...')
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
        print('Running DADA2 taxonomy classification ...')
        taxa = f'{TMP}/{base}_taxa.csv'
        dada2_taxa = dada2_taxonomy(asv, taxa, dada2_tax_dbs[taxreference])

        # join dada2 taxonomy results with ASV data
        output_data = pd.merge(
            asv_data,
            dada2_taxa,
            on='sequence',
            how='left'
        ).drop('index', axis=1)

    # output table with ASVs, abundance, and taxonomy
    output_data.to_csv(output, index=False)
    print(f'Analysis complete. Find results at {output}')

if __name__ == '__main__':
    main()
