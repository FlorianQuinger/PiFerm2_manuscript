#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=013_move_reads
#SBATCH --output=logs/013_move_reads.out
#SBATCH --error=logs/013_move_reads.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/

RR=01_raw_reads
UNTAR=01_raw_untar

mkdir -p $RR

for i in $UNTAR/*;do
 for j in $i/01.RawData/*;do
  echo $j
  mv $j $RR/
 done
done
