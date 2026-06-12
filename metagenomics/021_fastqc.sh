#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8g
#SBATCH --job-name=021_fastqc
#SBATCH --output=logs/021_fastqc.out
#SBATCH --error=logs/021_fastqc.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate FastQC_v0.12.1

#create directory
mkdir -p 02_report/02_fastqc/

#environment variables
RR=00_data/01_raw_reads
OUT=02_report/02_fastqc

#loop over files and create fastqc reports
for inp1 in $RR/*_1.fq.gz; do
 inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
 fastqc -t 2 -o $OUT $inp1 $inp2
 echo "processed $inp1 $inp2"
done

#create directory for multiqc
mkdir -p 02_report/02_multiqc/

multiqc $OUT -o 02_report/02_multiqc/



