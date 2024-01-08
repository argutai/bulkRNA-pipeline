#!/bin/bash -l
#SBATCH --job-name=featureCounts
#SBATCH --mem=80G
#SBATCH --output=logs/featureCounts/%j.featureCounts.out

allfiles=$(find . -name *.bam)
echo $allfiles

singularity run -B $(pwd) \
featurecount.sif featureCounts $allfiles \
-a gencode.v38.annotation.gtf \
-o results/featureCounts/RNA_seq_featureCounts.txt
