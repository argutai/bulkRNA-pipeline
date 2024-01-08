#!/bin/bash -l
#SBATCH --job-name=qualimap
#SBATCH --mem=80G
#SBATCH --output=logs/qualimap/%j.qualimap.out

SAMPLE=$1

#### QUALIMAP
singularity run -B $(pwd) \
qualimap.sif qualimap rnaseq \
-bam clean-S3/$SAMPLE/${SAMPLE}Aligned.out.bam \
-gtf gencode.v38.annotation.gtf \
-outdir results-S3/qualimap/$SAMPLE \
-outformat HTML \
--java-mem-size=30G

#### PICARD
singularity run -B $(pwd) \
picard.sif java -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
I=clean-S3/$SAMPLE/${SAMPLE}Aligned.out.bam \
O=results-S3/picard/$SAMPLE \
REF_FLAT=refFlat.txt \
STRAND=NONE
