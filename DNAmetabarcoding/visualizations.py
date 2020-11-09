#!/usr/bin/env python
import os
import sys
import glob
import pandas as pd
import subprocess

def run_phyloseq(ASVfile, taxafile, taxalevel):
    subprocess.check_call(['Rscript', 'phyloseq.R', ASVfile, taxafile, taxalevel])
    
def mergefiles(directory, taxalevel):
    files = glob.glob(f'{directory}/*.csv')

    # abundance data by sequence for each sample
    abundance = (
        pd.concat(
            [
                pd.read_csv(f)
                .set_index('sequence')
                [['abundance']]
                .rename(columns={'abundance': f.rsplit('.', 1)[0]})
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
    abundance.to_csv('mergedASVdata.csv')
    
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
    )
    taxa.to_csv('mergedtaxadata.csv')
    
    # check the lengths of the abundance and taxa data are the same
    # if not, this indicates that the same sequence was identified as
    # belonging to different taxa, and some of the data should be updated
    assert len(abundance) == len(taxa)

    run_phyloseq('mergedASVdata.csv', 'mergedtaxadata.csv', taxalevel)
    

if __name__ == "__main__":
    location = sys.argv[1]
    taxalevel = sys.argv[2]
    mergefiles(location, taxalevel)

