#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=011_unzip
#SBATCH --output=logs/011_unzip.out
#SBATCH --error=logs/011_unzip.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

#create directory
mkdir -p 00_data/01_raw_untar/

#environment variables
TAR=00_data/01_raw
UNTAR=00_data/01_raw_untar

#untar all .tar files
for i in $TAR/*.tar;do 
 echo $i
 tar -xvf $i -C $UNTAR
done





