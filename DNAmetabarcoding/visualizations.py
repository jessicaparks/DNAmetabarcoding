#!/usr/bin/env python
import os
import sys
import glob
import pandas as pd


def mergefiles(directory):
    directory = directory + r"\*.csv"
    #print(directory)

    abundancefilelist = []
    #print(type(abundancefilelist))

    #print(glob.glob(directory))
    #df = pd.DataFrame()
    for counter, files in enumerate(glob.glob(directory)):
        #print(counter)
        if counter == 0:
            taxafile = pd.read_csv(files, index_col = 0, usecols = [0,2,3,4,5,6,7,8])
        else:
            taxafile = taxafile.merge(pd.read_csv(files, index_col = 0, usecols = [0,2,3,4,5,6,7,8]), on = "sequence", how = 'outer')
            #missing_x = taxafile[taxafile['kingdom_x'].empty].index.tolist()
            #if taxafile.kingdom_x.to_string(index=False) == "" and taxafile.kingdom_y.to_string() 
            #print(missing_x)

            #I was attempting to do some type of filtering here in the comments to get it back to one set of columns but including the new data for each merge
        abundancefilelist.append(pd.read_csv(files, index_col = 0, usecols = [0,1]))
    #print(taxafile)
    abundancedataframe = pd.concat(abundancefilelist, axis = 1)
    abundancedataframe.to_csv("mergedASVdata.csv")
    taxafile.to_csv("mergedtaxadata.csv")

if __name__ == "__main__":
    location = sys.argv[1]
    mergefiles(location)

