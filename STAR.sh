#!/bin/bash -l
#SBATCH --job-name=STAR
#SBATCH --output=logs/clean/%j.STAR.out
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=80G
#SBATCH --time=0-48:00

SAMPLE=$1
mkdir clean-S3/$SAMPLE	

reference=index/genome/STAR
STAR-2.7.11a/bin/Linux_x86_64/STAR \
	--outSAMtype BAM Unsorted \
	--runThreadN 30 \
	--readFilesCommand zcat \
	--genomeDir $reference \
	--readFilesIn path/to/data/${SAMPLE}/${SAMPLE}_1.fq.gz path/to/data/${SAMPLE}/${SAMPLE}_2.fq.gz \
	--outFileNamePrefix clean-S3/${SAMPLE}/${SAMPLE}
