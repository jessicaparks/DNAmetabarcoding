#!/bin/tcsh

# activate the conda environment
conda activate /usr/local/usrapps/trnL_blast/conda/envs/dnametabarcoding

# set the input list file,  output directory, the primer fasta file, and the NCBI entrez key file
set filelist = /share/trnL_blast/data/ITS2_list.txt
set outdir = /share/trnL_blast/data/ITS2_results
set primers = /share/trnL_blast/data/ITS2_primers.fasta
mkdir -p $outdir

# submit a job to the cluster for each input file
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
    "python app.py -i ${file} -o ${outdir}/${filename}.csv --primers ${primers} --taxmethod BLAST --threads 8"
end

# deactivate the conda environment
conda deactivate
