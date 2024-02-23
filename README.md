# BulkRNA Pipeline and Analysis

Pull images (use srun if on CREATE HPC):
```
singularity pull fastqc.sif docker://alexcoppe/fastqc
singularity pull picard.sif docker://broadinstitute/picard
singularity pull featurecount.sif docker://alexgilgal/featurecount
```
Install [STAR](https://github.com/alexdobin/STAR), add to path

### Fast QC
```
for FILE_PATH in $(find path/to/data -name *.fq.gz); do
	echo $FILE_PATH
	sbatch -p cpu fastqc.sh $FILE_PATH
done
```

### STAR alignment
```
for SAMPLE in $(ls path/to/data); do
	echo $SAMPLE
	sbatch -p cpu STAR.sh $SAMPLE
done
```

### Qualimap & Picard
```
for SAMPLE in $(ls path/to/data); do
	echo $SAMPLE
	sbatch -p cpu picard.sh $SAMPLE
done
```

### Feature maps

```
sbatch -p cpu featurecounts.sh
```
