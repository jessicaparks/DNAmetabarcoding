# DNAmetabarcoding


### Setting up your environment on Henry2
The DNAmetabarcoding program uses a conda environment. To enable conda access on the HPC, you should run the following two commands and then log out and log back in. This will only need to be done once. For more details on conda setup on use on Henry2, see https://projects.ncsu.edu/hpc/Software/Apps.php?app=Conda. You may need to complete additional steps from these instructions if you want to create or edit the Conda environment.  
```bash
module load conda
conda init tcsh
```

### The DNAmetabarcoding Conda environment


### Running the DNAmetabarcoding program


### Updating the nt BLAST database
By default, this program uses the nt BLAST database hosted by Henry2. Documentation for BLAST databases available on Henry2 is found at https://projects.ncsu.edu/hpc/Software/Apps.php?app=BLAST. The nt database is found at `/gpfs_partners/databases/ncbi/blast/nt/nt`.  
Users can optionally install and use their own copy of the nt database. To download the newest version of the nt BLAST database from NCBI, move to the `/share/trnL_blast/ncbi-reference/` directory and run the following commands.  
```bash
cd /share/trnL_blast/ncbi-reference
conda activate dnametabarcoding
update_blastdb.pl --decompress nt
```
The newest versions of the NCBI BLAST databases can be viewed at https://ftp.ncbi.nlm.nih.gov/blast/db/.
