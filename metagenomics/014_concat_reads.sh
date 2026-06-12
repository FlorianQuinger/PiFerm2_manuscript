#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=014_concat_reads
#SBATCH --output=logs/014_concat_reads.out
#SBATCH --error=logs/014_concat_reads.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/01_raw_reads

#concat reads with more than one file per forward/backward, rename, and move
for i in *;do
 echo $i
 inp1=$i/*_1.fq.gz
 inp2=$i/*_2.fq.gz
 echo "concatenating $inp1 $inp2"
 out1=${i}_1.fq.gz
 out2=${i}_2.fq.gz
 echo "moved $out1 $out2"
 cat $inp1 > $out1
 cat $inp2 > $out2
 echo "removing $i"
 rm -r $i
done
