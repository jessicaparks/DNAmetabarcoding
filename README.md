# DNAmetabarcoding

### Setting up your environment on Henry2
The DNAmetabarcoding program uses a conda environment. To enable conda access on the HPC, you should run the following two commands and then log out and log back in. This will only need to be done once. For more details on conda setup on use on Henry2, see https://projects.ncsu.edu/hpc/Software/Apps.php?app=Conda.  
```bash
module load conda
conda init tcsh
```
You will need to complete this additional step (also detailed in the instructions at https://projects.ncsu.edu/hpc/Software/Apps.php?app=Conda) if you want to create or edit the Conda environment. Create a file in your home directory named `.condarc`, with the following contents:  
```bash
channels:
  - bioconda
  - conda-forge
  - defaults
envs_dirs:
  - /usr/local/usrapps/trnL_blast/conda/envs
pkgs_dirs:
  - /usr/local/usrapps/trnL_blast/conda/pkgs
report_errors: false
```
This will direct conda to store environments and packages at `/usr/local/usrapps/trnL_blast/conda/`. If you do not have this `.condarc` file, conda will attempt to store these files in your home directory, which has very limited space on the HPC. This will result in an error about running out of space.  

### The DNAmetabarcoding Conda environment
To list the existing conda environments, run this command:  
```bash
conda env list
```
To create the conda environment, you can run the following command:  
(This should either be run from the location of the `environment.yaml` file or use the path to that file.)  
```bash
conda env create -f environment.yaml
```
To update the conda environment, you can run this command:  
```bash
conda env update -f environment.yaml
```

### Running the DNAmetabarcoding program




### Running the visualizations




### Updating the NCBI Taxize database




### Updating the nt BLAST database
By default, this program uses the nt BLAST database hosted by Henry2. Documentation for BLAST databases available on Henry2 is found at https://projects.ncsu.edu/hpc/Software/Apps.php?app=BLAST. The nt database is found at `/gpfs_partners/databases/ncbi/blast/nt/nt`.  
Users can optionally install and use their own copy of the nt database. To download the newest version of the nt BLAST database from NCBI, move to the `/usr/local/usrapps/trnL_blast/ncbi` directory and run the following commands.  
```bash
cd /usr/local/usrapps/trnL_blast/ncbi
conda activate /usr/local/usrapps/trnL_blast/conda/envs/dnametabarcoding
update_blastdb.pl --decompress nt
```
The newest versions of the NCBI BLAST databases can be viewed at https://ftp.ncbi.nlm.nih.gov/blast/db/.

### Downloading data from the HPC cluster
