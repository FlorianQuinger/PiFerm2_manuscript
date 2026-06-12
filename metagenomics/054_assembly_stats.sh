#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=64g
#SBATCH --job-name=054_assembly_stats
#SBATCH --output=logs/054_assembly_stats.out
#SBATCH --error=logs/054_assembly_stats.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate bbmap_v39.17

#environment variables
BR=00_data/04_bins_renamed
OUT=02_report/05_assembly_stats

#create directory
mkdir -p $OUT

# create assembly stats

statswrapper.sh $BR/*.fna > $OUT/stats.txt



