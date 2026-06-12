#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --mem=120g
#SBATCH --job-name=031_megahit
#SBATCH --output=logs/031_megahit.out
#SBATCH --error=logs/031_megahit.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate megahit_v1.2.9

#create directory
mkdir -p 00_data/03_contigs_megahit # megahit requires to create the output directory by itself

#environment variables
RC=00_data/02_reads_clean
MC=00_data/03_contigs_megahit

# loop megahit over all samples
for inp1 in $RC/*_1.fq.gz; do
	inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
	# specify directory for each sample
	sampleid=$(basename "${inp1/_1.fq.gz/}")
	if [ -d "$MC/$sampleid" ]; then #check if directory already exists, if continuing script
		echo "$sampleid allready processed"
	else
		megahit -1 $inp1 -2 $inp2 -o $MC/$sampleid --k-min 27 --k-max 97 --k-step 10
		echo "processed sample $sampleid"
	fi
done
