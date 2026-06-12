#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=024_kraken_unzip
#SBATCH --output=logs/024_kraken_unzip.out
#SBATCH --error=logs/024_kraken_unzip.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/00_databases

# download database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_gtdb_genome_reps_20241109.tar.gz

#create directory
mkdir -p kraken2/

#environment variables
TAR=k2_gtdb_genome_reps_20241109.tar.gz
UNTAR=kraken2

#untar all .tar file
tar -xvf $TAR -C $UNTAR
