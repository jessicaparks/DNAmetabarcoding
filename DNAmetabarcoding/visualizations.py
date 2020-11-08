#!/usr/bin/env python
import os
import sys
import glob
import pandas as pd
import subprocess

def run_phyloseq(ASVfile, taxafile, taxalevel):
    subprocess.check_call(['Rscript', 'phyloseq.R', ASVfile, taxafile, taxalevel])
    
def mergefiles(directory, taxalevel):
    #directory = directory + r"\*.csv"
    files = glob.glob(f'{directory}/*.csv')
    #print(directory)

    #abundancefilelist = []
    #print(type(abundancefilelist))

    #print(glob.glob(directory))
    #df = pd.DataFrame()
    
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
    #for counter, files in enumerate(glob.glob(directory)):
        #print(counter)
    #    if counter == 0:
    #        taxafile = pd.read_csv(files, index_col = 0, usecols = [0,2,3,4,5,6,7,8])
    #    else:
    #        taxafile = taxafile.merge(pd.read_csv(files, index_col = 0, usecols = [0,2,3,4,5,6,7,8]), on = "sequence", how = 'outer')
            #missing_x = taxafile[taxafile['kingdom_x'].empty].index.tolist()
            #if taxafile.kingdom_x.to_string(index=False) == "" and taxafile.kingdom_y.to_string() 
            #print(missing_x)

            #I was attempting to do some type of filtering here in the comments to get it back to one set of columns but including the new data for each merge
    #    abundancefilelist.append(pd.read_csv(files, index_col = 0, usecols = [0,1]))
    #print(taxafile)
    #abundancedataframe = pd.concat(abundancefilelist, axis = 1)
    #abundancedataframe.to_csv("mergedASVdata.csv")
    #taxafile.to_csv("mergedtaxadata.csv")
    

if __name__ == "__main__":
    location = sys.argv[1]
    taxalevel = sys.argv[2]
    mergefiles(location, taxalevel)

