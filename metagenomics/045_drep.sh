#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=045_drep
#SBATCH --output=logs/045_drep.out
#SBATCH --error=logs/045_drep.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate drep_v3.5.0

#create directory
mkdir -p 00_data/04_bins_drep
mkdir -p 00_data/04_bins_das_tool/all_bins
mkdir -p 02_report/04_bins_drep

#environment variables
BDT=00_data/04_bins_das_tool
BDR=00_data/04_bins_drep
BDRout=02_report/04_bins_drep

# copy all bins to new directory and name uniquely

for FOLDER in $BDT/*; do
        sampleid=$(basename $FOLDER)
	for FILE in $FOLDER/das_tool_DASTool_bins/*.fa; do
		filename=$(basename $FILE)
	  	cp $FILE $BDT/all_bins/${sampleid}_${filename}
	done
done

# drep
dRep dereplicate $BDR -g $BDT/all_bins/*.fa -sa 0.98 -p 32

# move drep plots to reports folder

cp -r $BDR/figures $BDRout
