#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=065_pig_gut_prokka
#SBATCH --output=logs/065_pig_gut_prokka.out
#SBATCH --error=logs/065_pig_gut_prokka.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate prokka_v1.14.6

#environment variables
PGF=00_data/00_databases/pig_gut/fna_files
PGP=00_data/06_pig_gut_prokka

#create directory
mkdir -p $PGP

# copy all bins to new directory

for FILE in $PGF/*; do
	MAG=$(basename "$FILE")
	MAG=${MAG/.fna/}
	prokka --outdir $PGP --prefix $MAG --locustag $MAG --cpus 32 --force $FILE 
	echo "processed $MAG"
done

