#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=062_prokka
#SBATCH --output=logs/062_prokka.out
#SBATCH --error=logs/062_prokka.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate prokka_v1.14.6

#environment variables
BR=00_data/04_bins_renamed
BP=00_data/06_bins_prokka

#create directory
mkdir -p $BP

# copy all bins to new directory

for FILE in $BR/*; do
	MAG=$(basename "$FILE")
	MAG=${MAG/.fna/}
	prokka --outdir $BP --prefix $MAG --locustag $MAG --cpus 32 --force $FILE 
	echo "processed $MAG"
done

