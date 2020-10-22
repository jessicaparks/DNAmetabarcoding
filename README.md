# DNAmetabarcoding


### Setting up your environment on Henry2
The DNAmetabarcoding program uses a conda environment. To enable conda access on the HPC, you should run the following two commands and then log out and log back in. This will only need to be done once. For more details on conda setup on use on Henry2, see https://projects.ncsu.edu/hpc/Software/Apps.php?app=Conda.  
```bash
module load conda
conda init tcsh
```

### Updating the nt BLAST database
To download the newest version of the nt BLAST database from NCBI, move to the `/share/trnL_blast/ncbi-reference/` directory and run the following commands.  
```bash
cd /share/trnL_blast/ncbi-reference
conda activate dnametabarcoding
update_blastdb.pl --num_threads 1 --decompress nt
```
The newest versions of the NCBI BLAST databases can be found at https://ftp.ncbi.nlm.nih.gov/blast/db/.
