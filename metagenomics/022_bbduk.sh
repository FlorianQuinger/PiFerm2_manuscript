#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --time=36:00:00
#SBATCH --mem=100g
#SBATCH --job-name=022_bbduk
#SBATCH --output=logs/022_bbduk.out
#SBATCH --error=logs/022_bbduk.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate bbmap_v39.17

#create directory
mkdir -p 00_data/02_reads_trimmed
mkdir -p 02_report/02_fastqc_trimmed
mkdir -p 00_data/02_reads_filtered
mkdir -p 02_report/02_fastqc_filtered
mkdir -p 00_data/02_reads_qtrimmed
mkdir -p 02_report/02_fastqc_qtrimmed
mkdir -p 02_report/02_multiqc_bbduk

#environment variables
RR=00_data/01_raw_reads
RT=00_data/02_reads_trimmed
RF=00_data/02_reads_filtered
RQ=00_data/02_reads_qtrimmed
OUTtrim=02_report/02_fastqc_trimmed
OUTfilt=02_report/02_fastqc_filtered
OUTqtrim=02_report/02_fastqc_qtrimmed
BBMAP_RESOURCES=/home/ho/ho_genetics/ho_yogau97/.conda/envs/bbmap_v39.17/share/bbmap/resources

#loop over files and create fastqc reports
for inp1 in $RR/*_1.fq.gz; do
#for inp1 in $RR/AP2_101_1.fq.gz; do # only one sample
 inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
 outp1=${inp1/${RR}/${RT}}
 outp2=${inp2/${RR}/${RT}}
 
 #Trimming
 conda activate bbmap_v39.17 # activate in case not activated (loop)
 bbduk.sh in1=$inp1 in2=$inp2 out1=$outp1 out2=$outp2 ref=$BBMAP_RESOURCES/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe
 conda activate FastQC_v0.12.1 # switch environment
 fastqc -t 2 -o $OUTtrim $outp1 $outp2
 echo "trimmed $inp1 $inp2"
 #rm $inp1 $inp2 # not removing 
 
 #Contaminant filtering
 conda activate bbmap_v39.17 # switch back environment
 inp1=$outp1
 inp2=$outp2
 outp1=${inp1/${RT}/${RF}}
 outp2=${inp2/${RT}/${RF}}
 bbduk.sh in1=$inp1 in2=$inp2 out1=$outp1 out2=$outp2 k=31 ref=$BBMAP_RESOURCES/sequencing_artifacts.fa.gz,$BBMAP_RESOURCES/phix174_ill.ref.fa.gz
 conda activate FastQC_v0.12.1 # switch environment
 fastqc -t 2 -o $OUTfilt $outp1 $outp2
 echo "filtered $inp1 $inp2"
 rm $inp1 $inp2

 #Quality trimming and filtering
 conda activate bbmap_v39.17 # switch back environment
 inp1=$outp1
 inp2=$outp2
 outp1=${inp1/${RF}/${RQ}}
 outp2=${inp2/${RF}/${RQ}}
 bbduk.sh in1=$inp1 in2=$inp2 out1=$outp1 out2=$outp2 minlength=51 qtrim=rl maq=10 trimq=10 maxns=0 #minlength 1/3 of read length
 conda activate FastQC_v0.12.1 # switch environment
 fastqc -t 2 -o $OUTqtrim $outp1 $outp2
 echo "quality trimmed $inp1 $inp2"
 rm $inp1 $inp2
done

#create multiqc reports for each level

multiqc 02_report/02_fastqc_trimmed/ -o 02_report/02_multiqc_bbduk/
multiqc 02_report/02_fastqc_filtered/ -o 02_report/02_multiqc_bbduk/
multiqc 02_report/02_fastqc_qtrimmed/ -o 02_report/02_multiqc_bbduk/



