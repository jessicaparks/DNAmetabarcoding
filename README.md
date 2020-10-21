# DNAmetabarcoding


### Updating the nt BLAST database
To download the newest version of the nt BLAST database from NCBI, move to the `/share/trnL_blast/ncbi-reference/` directory and run the following commands.  
```bash
cd /share/trnL_blast/ncbi-reference
conda activate dnametabarcoding
update_blastdb.pl --num_threads 1 --decompress nt
```
The newest versions of the NCBI BLAST databases can be found at https://ftp.ncbi.nlm.nih.gov/blast/db/.
