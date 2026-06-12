#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=00_checksum
#SBATCH --output=logs/00_checksum.out
#SBATCH --error=logs/00_checksum.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/01_raw

md5sum -c MD5\ \(1\).txt





