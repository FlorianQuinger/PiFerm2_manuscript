#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=033_coverm
#SBATCH --output=logs/033_coverm.out
#SBATCH --error=logs/033_coverm.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate coverm_v0.7.0

#create directory
mkdir -p 00_data/03_contigs_coverm
mkdir -p 00_data/03_contigs_bam

#environment variables
RC=00_data/02_reads_clean
FC=00_data/03_contigs_filtered
CC=00_data/03_contigs_coverm
CB=00_data/03_contigs_bam

# loop over all read files
for inp1 in $RC/*_1.fq.gz; do 
	inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
	#define input directory for contigs
	FOLDER=$FC/$(basename "${inp1/_1.fq.gz/}")
	#create output directories
	OUTFOLDERcoverm=$CC/$(basename $FOLDER)
	OUTFOLDERbam=$CB/$(basename $FOLDER)
	mkdir -p $OUTFOLDERcoverm
	mkdir -p $OUTFOLDERbam
	#CoverM
	coverm contig -m metabat -1 $inp1 -2 $inp2 -r $FOLDER/final.contigs.fa -o $OUTFOLDERcoverm/coverm_depth.txt --bam-file-cache-directory $OUTFOLDERbam/ -t 32 
	echo "processed $FOLDER"
done
