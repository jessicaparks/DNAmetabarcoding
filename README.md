# DNAmetabarcoding
Identification and taxonomic classification of exact amplicon sequence variants (ASVs) from DNA metabarcoding data.

## Table of Contents
* [Installation](#installation)
* [Quick Start](#quick-start)
  * [single sample](#analyze-a-single-sample)
  * [multiple samples](#analyze-multiple-samples)
  * [visualization](#create-an-abundance-plot)
* [Environment](#environment)
  * [Setting up your environment on Henry2](#setting-up-your-environment-on-henry2)
  * [The DNAmetabarcoding Conda environment](#the-dnametabarcoding-conda-environment)
* [The DNAmetabarcoding program](#the-dnametabarcoding-program)
* [Visualization](#visualization)
* [Databases](#databases)
  * [DADA2 taxonomy databases](#dada2-taxonomy-databases)
  * [NCBI nt BLAST database](#ncbi-nt-blast-database)
  * [Taxize NCBI database](#taxize-ncbi-database)
* [Data Storage and Retrieval](#data-storage-and-retrieval)
  * [Downloading data from the HPC cluster](#downloading-data-from-the-hpc-cluster)

## Installation
To install the DNAmetabarcoding program, use the following command. This path can also be copied from the `Code` dropdown menu.
```bash
git clone https://github.com/jessicaparks/DNAmetabarcoding.git
```
This will download the code to a directory named `DNAmetabarcoding`. For authentication, it will ask for your GitHub username and password.

## Quick Start
Each of the following commands should be run from within the DNAmetabarcoding sub-directory. Before running any of these commands, you'll need to do the first step under [Setting up your environment on Henry2](#setting-up-your-environment-on-henry2) to install Conda.
#### Analyze a single sample:
Run the following command with the arguments substituted to match your data.
```bash
conda activate /usr/local/usrapps/trnL_blast/conda/envs/dnametabarcoding
./app.py -i INPUT_FASTQ_FILE -o OUTPUT_CSV_FILE --primers PRIMER_FASTA_FILE --taxmethod BLAST
```

#### Analyze multiple samples:
Edit the `submit_job.sh` script and run `./submit_job.sh` to send the jobs to the cluster.

#### Create an abundance plot:
Run the following command with the correct arguments for your data. This requires that the csv files that you wish to combine and visualize are in a directory together.
```bash
conda activate /usr/local/usrapps/trnL_blast/conda/envs/dnametabarcoding
./visualizations.py -d INPUT_DIRECTORY -o OUTPUT_DIRECTORY -p OUTPUT_FILE_PREFIX -r TAXONOMIC_RANK -f TOP_TAXA_FILTER
```

## Environment
The DNAmetabarcoding program includes an environment that is managed with [Conda](https://docs.conda.io/en/latest/). Instructions for working with and editing this environment are in the sections below.

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

## The DNAmetabarcoding program




## Visualization



## Databases

### DADA2 taxonomy databases

### NCBI nt BLAST database
By default, this program uses the nt BLAST database hosted by Henry2. Documentation for BLAST databases available on Henry2 is found at https://projects.ncsu.edu/hpc/Software/Apps.php?app=BLAST. The nt database is found at `/gpfs_partners/databases/ncbi/blast/nt/nt`.  
Users can optionally install and use their own copy of the nt database. To download the newest version of the nt BLAST database from NCBI, move to the `/usr/local/usrapps/trnL_blast/ncbi` directory and run the following commands.  
```bash
cd /usr/local/usrapps/trnL_blast/ncbi
conda activate /usr/local/usrapps/trnL_blast/conda/envs/dnametabarcoding
update_blastdb.pl --decompress nt
```
The newest versions of the NCBI BLAST databases can be viewed at https://ftp.ncbi.nlm.nih.gov/blast/db/.

### Taxize NCBI database

## Data Storage and Retrieval

### Downloading data from the HPC cluster

