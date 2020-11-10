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
    """
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
	.sort_values(ascending=False)[:20]
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
	raise ValueError('ASVs are not uniquely identified from all samples. '
			 'Check that the samples were analyzed with the same taxonomy '
			 'method (BLAST or DADA2) and the same database.')

    plot_file = f'{outputdir}/{outputprefix}_abundance_plot.png'
    run_phyloseq(abundance_file, filt_taxa_file, rank, plot_file)


if __name__ == '__main__':
    main()
