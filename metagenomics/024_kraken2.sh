#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=336:00:00
#SBATCH --mem=800g
#SBATCH --job-name=024_kraken2
#SBATCH --output=logs/024_kraken2.out
#SBATCH --error=logs/024_kraken2.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate kraken2_v2.1.4

#create directory
mkdir -p 02_report/02_kraken2

#environment variables
RC=00_data/02_reads_clean
OUT=02_report/02_kraken2
DB=00_data/00_databases/kraken2

# run kraken2 per sample
for inp1 in $RC/*_1.fq.gz; do 
 inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
 # sampleid
 sampleid=$(basename "${inp1/_1.fq.gz/}")
 echo "start kraken2 for $sampleid"
 kraken2 --db $DB --threads 48 --report $OUT/"$sampleid".k2report --paired --minimum-hit-groups 3 --confidence 0.5 $inp1 $inp2 > $OUT/"$sampleid".kraken2

 echo "start bracken for $sampleid"
 bracken -d $DB -i $OUT/"$sampleid".k2report -r 100 -l S -t 10 -o $OUT/"$sampleid".bracken -w $OUT/"$sampleid".breport #read length set to approximately shortest read length observed
 
 echo "processed $inp1 $inp2"
done
