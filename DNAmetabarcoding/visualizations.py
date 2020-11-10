#!/usr/bin/env python

import click
import glob
import os
import pandas as pd
import subprocess

# set options for taxonomic ranks
ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def run_phyloseq(asvfile, taxafile, taxarank, plotfile):
    """Run phyloseq to create an abundance plot.
    The file will be output to the location specified by plotfile.
    
    Arguments:
    asvfile -- a csv file containing the ASVs and their abundance in each sample
    taxafile -- a csv file containing the ASVs and their taxonomy
    taxarank -- the taxonomic rank at which to analyze the data
    plotfile -- output path for the figure (png format)
    """
    subprocess.check_call(['Rscript', 'phyloseq.R', asvfile, taxafile, taxarank, plotfile])


@click.command()
@click.option('-d', '--directory', type=click.Path(exists=True, file_okay=False),
	      required=True, help='path for directory containing the input CSV files')
@click.option('-o', '--outputdir', type=click.Path(file_okay=False), required=True,
	      help='path for the output directory')
@click.option('-p', '--outputprefix', type=str, required=True,
	      help='prefix for the output files')
@click.option('-r', '--rank', type=click.Choice(ranks, case_sensitive=False),
	      required=True, help='taxonomic rank to be visualized')
@click.option('-f', '--filter', type=int, default=10, show_default=True,
	      help='number of top taxa to display in visualization')
def main(directory, outputdir, outputprefix, rank, filter):
    """Summarize the ASV abundance and taxonomy data and plot the abundance.
    The ASV abundance and taxonomy data from all of the CSV files in the input
    directory is merged by the ASV sequence.

    This will produce the following output files in the specified output directory:
    (1) a csv with all the taxonomy data, (2) a csv with all the abundance data,
    (3) a csv with the merged taxonomy and abundance data, (4) a csv with the
    filtered taxonomy data where the taxon of interest is replaced with "Other"
    for values below the cutoff rank abundance, and (5) a png-formated image file
    with the abundance plot.

    "Unknown" indicates taxonomy values that were not determined. "Other" indicates
    taxonomy values that are below the rank abundance cutoff, when considering
    cumulative abundance in all samples in the input directory.

    Note: This relies on the taxonomy determined for each of the samples being
    consistent for the same ASV sequence. If different methods (BLAST and DADA2)
    were used to determine the taxonomy, this likely will not be true and will
    raise an error.\f

    Arguments:
    directory -- path for directory containing the input CSV files
    outputdir -- path for the output directory
    outputprefix -- prefix for the output files
    rank -- taxonomic rank to be visualized
    filter -- number of top taxa to display in visualization
    """
    # read the list of CSV files from the input directory
    files = glob.glob(f'{directory}/*.csv')
    # create the output directory if it does not already exist
    outputdir = outputdir.rstrip('/')
    os.makedirs(outputdir, exist_ok=True)

    # abundance data by sequence for each sample
    abundance = (
        pd.concat(
            [
             	pd.read_csv(f)
                .set_index('sequence')
                [['abundance']]
                .rename(columns={'abundance': os.path.basename(f).rsplit('.', 1)[0]})
                for f in files
            ],
            axis=1,
            sort=False
        )
        .fillna(0)
        .astype(int)
        .reset_index()
        .rename(columns={'index': 'sequence'})
    )
    abundance_file = f'{outputdir}/{outputprefix}_ASV.csv'
    abundance.to_csv(abundance_file, index=False)
    print(f'Abundance data: {abundance_file}')

    # taxonomy classification data for each sequence
    taxa = (
        pd.concat(
            [
             	pd.read_csv(f)
                [['sequence', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']]
                for f in files
            ],
            axis=0,
            sort=False,
            ignore_index=True
        )
        .drop_duplicates()
        .reset_index(drop=True)
        .fillna('Unknown')
    )
    taxa_file = f'{outputdir}/{outputprefix}_taxa.csv'
    taxa.to_csv(taxa_file, index=False)
    print(f'Taxonomy data: {taxa_file}')

    # create a merged ASV and taxa file for all the samples
    all_data = pd.merge(
        abundance,
        taxa,
        on='sequence',
        how='inner'
    )
    all_data_file = f'{outputdir}/{outputprefix}_all.csv'
    all_data.to_csv(all_data_file, index=False)
    print(f'Merged ASV, abundance, and taxonomy data: {all_data_file}')

    # filter to top most abundant taxa
    print(f'Filtering {rank} to the top {filter} most abundant ...')
    count_remove = len(all_data[all_data[rank]!='Unknown'][rank].unique()) - filter
    print(f'{count_remove} {rank} values will be replaced by "Other".')
    keep = (
        all_data[all_data[rank]!='Unknown']
        .groupby(rank).sum()
        .sum(axis=1)
        .sort_values(ascending=False)[:filter]
        .index.tolist() + ['Unknown']
    )
    taxa[rank] = taxa[rank].apply(lambda x: x if x in keep else 'Other')
    
    # write taxa files for filtered taxa
    filt_taxa_file = f'{outputdir}/{outputprefix}_filtered_taxa.csv'
    taxa.to_csv(filt_taxa_file, index=False)
    print(f'Filtered taxonomy data: {filt_taxa_file}')

    # check the lengths of the abundance and taxa data are the same
    # if not, this indicates that the same sequence was identified as
    # belonging to different taxa, and some of the data should be updated
    if len(abundance) != len(taxa):
        raise ValueError(
            'ASVs are not uniquely identified from all samples. '
            'Check that the samples were analyzed with the same taxonomy '
            'method (BLAST or DADA2) and the same database.'
    )

    plot_file = f'{outputdir}/{outputprefix}_abundance_plot.png'
    print('Creating plot ...')
    run_phyloseq(abundance_file, filt_taxa_file, rank, plot_file)
    print(f'Abundance plot: {plot_file}')


if __name__ == '__main__':
    main()
