#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=044_das_tool
#SBATCH --output=logs/044_das_tool.out
#SBATCH --error=logs/044_das_tool.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate das_tool_v1.1.7

#create directory
mkdir -p 00_data/04_bins_das_tool

#environment variables
DT=$HOME/miniconda3/envs/das_tool_v1.1.7/bin
FC=00_data/03_contigs_filtered
BM=00_data/04_bins_metabat2
BS=00_data/04_bins_semibin2
BC=00_data/04_bins_comebin_cuda
BDT=00_data/04_bins_das_tool

# loop over all folders created by assembly

for FOLDER in $FC/*; do
	sampleid=$(basename $FOLDER)
	# create output directory
	OUTFOLDER=$BDT/$sampleid
	mkdir -p $OUTFOLDER
	#Preprocessing: create Contigs2Bin file
	Fasta_to_Contig2Bin.sh -e fa -i $BM/$sampleid > $OUTFOLDER/contigs2bin_metabat2.tsv
	Fasta_to_Contig2Bin.sh -e fa -i $BS/$sampleid/output_bins > $OUTFOLDER/contigs2bin_semibin2.tsv
	Fasta_to_Contig2Bin.sh -e fa -i $BC/$sampleid/comebin_res/comebin_res_bins > $OUTFOLDER/contigs2bin_comebin.tsv
	# DAS_Tool command
	DAS_Tool -i $OUTFOLDER/contigs2bin_metabat2.tsv,$OUTFOLDER/contigs2bin_semibin2.tsv,$OUTFOLDER/contigs2bin_comebin.tsv -c $FOLDER/final.contigs.fa -o $OUTFOLDER/das_tool -l metabat2,semibin2,comebin -t 32 --write_bin_evals --write_bins
	echo "processed sample $sampleid"
done


