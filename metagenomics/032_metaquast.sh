#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=100g
#SBATCH --job-name=032_metaquast
#SBATCH --output=logs/032_metaquast.out
#SBATCH --error=logs/032_metaquast.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate quast_v5.3.0

#create directory
mkdir -p 02_report/03_contigs_quast

#environment variables
FC=00_data/03_contigs_filtered
OUTmegahit=02_report/03_contigs_quast

#Metaquast for megahit
for FOLDER in $FC/*; do
	OUTFOLDER=$OUTmegahit/$(basename $FOLDER)
	mkdir -p $OUTFOLDER
	metaquast.py $FOLDER/final.contigs.fa -o $OUTFOLDER --max-ref-number 0
	echo "saved in $OUTFOLDER"
done
