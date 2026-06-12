#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:a30:1
#SBATCH --time=168:00:00
#SBATCH --mem=120g
#SBATCH --job-name=043_comebin_cuda
#SBATCH --output=logs/043_comebin_cuda.out
#SBATCH --error=logs/043_comebin_cuda.err

module load devel/miniforge
conda activate comebin_v1.0.4

cd ~/.conda/envs/comebin_v1.0.4/bin/COMEBin # working from comebin directory!

#environment variables
WORK=/pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics
RC=$WORK/00_data/02_reads_clean
FC=$WORK/00_data/03_contigs_filtered
BC=$WORK/00_data/04_bins_comebin_cuda

#create directory
mkdir -p $BC
mkdir -p $BC/coverage

for inp1 in $RC/*_1.fq.gz; do 
	inp2=${inp1/_1.fq.gz/_2.fq.gz} #replace _1 with _2
	# specify directory for each sample
	sampleid=$(basename "${inp1/_1.fq.gz/}")
	
	if [ -d "$BC/$sampleid" ]; then # check if folder already exists
		echo "$sampleid allready processed"
	else
		# create sample individual folder for coverage files
		mkdir -p $BC/coverage/$sampleid
		#create mapping files
		scripts/gen_cov_file.sh -a $FC/$sampleid/final.contigs.fa -o $BC/coverage/$sampleid -f "_1.fq.gz" -r "_2.fq.gz" -t 32 -l 2500 $RC/$sampleid*
		# create sample individual folders for bin output
		mkdir -p $BC/$sampleid
		# run comebin
		run_comebin.sh -a $FC/$sampleid/final.contigs.fa -o $BC/$sampleid -p $BC/coverage/$sampleid/work_files -t 32
		echo "processed sample $sampleid"
	fi
done
