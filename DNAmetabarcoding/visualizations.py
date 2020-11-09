#!/usr/bin/env python

import click
import os
import glob
import pandas as pd
import subprocess

# set options for taxonomic ranks
ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def run_phyloseq(ASVfile, taxafile, taxalevel):
    subprocess.check_call(['Rscript', 'phyloseq.R', ASVfile, taxafile, taxalevel])


@click.command()
@click.option('-d', '--directory', required=True)
@click.option('-o', '--output', required=True)
@click.option('-r', '--rank', type=click.Choice(ranks, case_sensitive=False), required=True)
def main(directory, output, rank):
    files = glob.glob(f'{directory}/*.csv')
    os.makedirs(output, exist_ok=True)

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
    abundance_file = f'{output}/mergedASVdata.csv'
    abundance.to_csv(abundance_file, index=False)

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
    taxa_file = f'{output}/mergedtaxadata.csv'
    taxa.to_csv(taxa_file, index=False)

    # check the lengths of the abundance and taxa data are the same
    # if not, this indicates that the same sequence was identified as
    # belonging to different taxa, and some of the data should be updated
    assert len(abundance) == len(taxa)

    run_phyloseq(abundance_file, taxa_file, rank)

if __name__ == '__main__':
    main()

