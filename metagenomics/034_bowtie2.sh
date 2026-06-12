#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=034_bowtie2
#SBATCH --output=logs/034_bowtie2.out
#SBATCH --error=logs/034_bowtie2.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

#to generate input files for Semibin2


module load devel/miniforge
#Mapping
conda activate bowtie2_v2.5.4

#environment variables
RC=00_data/02_reads_clean
FC=00_data/03_contigs_filtered
BT=00_data/03_contigs_bowtie2

#create directory
mkdir -p $BT
mkdir -p $BT/index
mkdir -p $BT/sam
mkdir -p $BT/bam
mkdir -p $BT/bam_filtered
mkdir -p $BT/bam_sorted


# loop over all read files
for inp1 in $RC/*_1.fq.gz; do
	inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
	sampleid=$(basename "${inp1/_1.fq.gz/}")
	#create output folder for index
	OUTindex=$BT/index/$sampleid
	bowtie2-build -f $FC/$sampleid/final.contigs.fa $BT/index/$sampleid
	#create output names for .sam
	OUTsam=$BT/sam/"$sampleid".sam
	#align with bowtie2
	bowtie2 -q --fr -x $BT/index/$sampleid -1 $inp1 -2 $inp2 -S $OUTsam -p 32
	# convert to bamfile
	OUTbam=$BT/bam/"$sampleid".bam
	samtools view -h -b -S $OUTsam -o $OUTbam -@ 32
	#filter alignments
	OUTfilter=$BT/bam_filtered/"$sampleid".filtered.bam
	samtools view -b -F 4 $OUTbam -o $OUTfilter -@ 32
	#sort bam file
	OUTsort=$BT/bam_sorted/"$sampleid".sorted.bam
	samtools sort -m 4000000000 $OUTfilter -o $OUTsort -@ 32

	echo "processed $sampleid"
done
