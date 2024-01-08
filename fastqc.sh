#!/bin/bash -l
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#BATCH --mem=20G
#SBATCH --output=logs/qc/%j.fastqc.out

FILE_PATH=$1
FILE=$FILE_PATH | sed 's:.*/::'
SAMPLE=$(echo $FILE | sed 's/.fq.gz//')

mkdir -p results/qc/$SAMPLE

singularity run -B /scratch/prj/cb_hormad1/ligand-bulkRNA \
fastqc.sif \
-o results/qc/$SAMPLE \
$FILE_PATH
