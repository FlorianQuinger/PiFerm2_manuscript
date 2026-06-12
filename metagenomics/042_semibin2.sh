#!/bin/bash

#SBATCH --partition=compute #not gpu since gpu was not able to be detected
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=64g
#SBATCH --job-name=042_semibin2
#SBATCH --output=logs/042_semibin2.out
#SBATCH --error=logs/042_semibin2.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate semibin2_v2.1.0

#create directory
mkdir -p 00_data/04_bins_semibin2

#environment variables
FC=00_data/03_contigs_filtered
BT=00_data/03_contigs_bowtie2
BS=00_data/04_bins_semibin2

# loop over read files to do single sample binning with pretrained model

for FOLDER in $FC/*; do
	# create output directory
	OUTFOLDER=$BS/$(basename $FOLDER)
	mkdir -p $OUTFOLDER
	# specify directory for each sample
	sampleid=$(basename $FOLDER)
	#semibin2
	SemiBin2 single_easy_bin -i $FOLDER/final.contigs.fa -b $BT/bam_sorted/${sampleid}.sorted.bam -o $OUTFOLDER --environment pig_gut --compression none --min-len 2500

	echo "processed sample $sampleid"
done
