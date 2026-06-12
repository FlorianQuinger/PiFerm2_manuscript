#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=120g
#SBATCH --job-name=041_metabat2
#SBATCH --output=logs/041_metabat2.out
#SBATCH --error=logs/041_metabat2.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate metabat2_v2.17

#create directory
mkdir -p 00_data/04_bins_metabat2

#environment variables
FC=00_data/03_contigs_filtered
CC=00_data/03_contigs_coverm
BM=00_data/04_bins_metabat2

# loop over contig folders
for FOLDER in $FC/*; do
	#path to coverm file
	COVERM=$CC/$(basename $FOLDER)
	OUTFOLDER=$BM/$(basename $FOLDER)
	mkdir -p $OUTFOLDER
	#Run metabat2
	metabat2 -i $FOLDER/final.contigs.fa -a $COVERM/coverm_depth.txt -o $OUTFOLDER/metabat2 -m 2500
	echo "processed $FOLDER"
done
