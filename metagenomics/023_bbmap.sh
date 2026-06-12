#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --time=128:00:00
#SBATCH --mem=300g
#SBATCH --job-name=023_bbmap
#SBATCH --output=logs/023_bbmap.out
#SBATCH --error=logs/023_bbmap.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate bbmap_v39.17

#create directory
mkdir -p 00_data/02_reads_contaminated
mkdir -p 00_data/02_reads_clean
mkdir -p 02_report/02_fastqc_clean

#environment variables
RQ=00_data/02_reads_qtrimmed
RC=00_data/02_reads_clean
RM=00_data/02_reads_contaminated
OUTclean=02_report/02_fastqc_clean
REF=00_data/00_databases/reference_masked.fa.gz

#loop over files and map to masked pig genome
for inp1 in $RQ/*_1.fq.gz; do
#for inp1 in $RQ/AP2_101_1.fq.gz; do # only one sample
 inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
 outpc1=${inp1/${RQ}/${RC}}
 outpc2=${inp2/${RQ}/${RC}}
 outpm1=${inp1/${RQ}/${RM}}
 outpm2=${inp2/${RQ}/${RM}}
 
 # map to pig reference genome
 conda activate bbmap_v39.17
 # parameters from bbmap/RemoveHuman
 bbmap.sh in1=$inp1 in2=$inp2 outm1=$outpm1 outm2=$outpm2 outu1=$outpc1 outu2=$outpc2 ref=$REF covstats=$OUTclean/covstats.txt minid=0.9 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 maxsites=1
 conda activate FastQC_v0.12.1 # change environment
 fastqc -t 2 -o $OUTclean $outpc1 $outpc2
 echo "mapped $inp1 $inp2"
 #rm $inp1 $inp2
done

#Multiqc report
multiqc 02_report/02_fastqc_clean -o 02_report/02_multiqc_bbduk/



