#!/bin/tcsh

# activate the conda environment
conda activate /usr/local/usrapps/trnL_blast/conda/envs/dnametabarcoding

# set the input list file, the primer fasta file, the output directory, and sequence cutoff size
set filelist = /share/trnL_blast/data/ITS2_list.txt
set primers = /share/trnL_blast/data/ITS2_primers.fasta
set outdir = /share/trnL_blast/data/ITS2_results
set cutoff = 50
# create the output directory
mkdir -p $outdir

# submit a job to the cluster for each input file
# the job name and outputs are named based on the filename (with .fastq.gz removed)
# the app.py command can be modified to run either BLAST or DADA2 taxonomy classification
foreach file ( `cat $filelist` )
  set filename = `basename $file .fastq.gz`
  bsub \
    -J $filename \
    -W 120 \
    -n 8 \
    -x \
    -R "span[hosts=1]" \
    -o $outdir/$filename.out.%J \
    -e $outdir/$filename.err.%J \
    "python app.py -i ${file} -o ${outdir}/${filename}.csv --primers ${primers} --taxmethod BLAST --threads 8 --cutoff ${cutoff}"
end

# deactivate the conda environment
conda deactivate
