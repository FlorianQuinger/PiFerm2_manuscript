#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=012_checksum
#SBATCH --output=logs/012_checksum.out
#SBATCH --error=logs/012_checksum.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/01_raw_untar

for i in *;do
 echo $i
 cd $i
 md5sum -c MD5.txt
 cd ..
done 




