#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=100g
#SBATCH --job-name=0312_filter_contigs
#SBATCH --output=logs/0312_filter_contigs.out
#SBATCH --error=logs/0312_filter_contigs.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate bbmap_v39.17

#create directory
mkdir -p 00_data/03_contigs_filtered

#environment variables
MC=00_data/03_contigs_megahit
FC=00_data/03_contigs_filtered

#loop over megahit folders and filter contigs for length 2500
for FOLDER in $MC/*; do
	OUTFOLDER=$FC/$(basename $FOLDER)
	mkdir -p $OUTFOLDER
	reformat.sh in=$FOLDER/final.contigs.fa out=$OUTFOLDER/final.contigs.fa minlength=2500
	echo "filtered $FOLDER"
done
