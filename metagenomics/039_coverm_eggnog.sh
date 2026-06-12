#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=039_coverm_eggnog
#SBATCH --output=logs/039_coverm_eggnog.out
#SBATCH --error=logs/039_coverm_eggnog.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate coverm_v0.7.0

#create directory
mkdir -p 02_report/03_contigs_coverm_eggnog

#environment variables
RC=00_data/02_reads_clean
EN=02_report/03_eggnog
CC=02_report/03_contigs_coverm_eggnog

# loop over all read files
for inp1 in $RC/*_1.fq.gz; do 
	inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
	#define input file for eggnog geneprediction
	FILE=$EN/$(basename "${inp1/_1.fq.gz/.emapper.genepred.fasta}") # replace with filename with geneprediction by eggnog
	sampleid=$(basename "${FILE/.emapper.genepred.fasta/}")
	#create genome definition file from geneprediction file
	awk -F' ' '/^>/ {print substr($1, 2) "\t" substr($1, 2)}' "$FILE" > $TMPDIR/${sampleid}.txt
	#CoverM
	coverm genome -m relative_abundance count -1 $inp1 -2 $inp2 -r $FILE --min-covered-fraction 0 --genome-definition $TMPDIR/${sampleid}.txt -o $CC/${sampleid}.txt -t 32 #according to https://github.com/wwood/CoverM/issues/211 
	echo "processed $FILE" 
done
