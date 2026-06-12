#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=061_db_drep
#SBATCH --output=logs/061_db_drep.out
#SBATCH --error=logs/061_db_drep.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate drep_v3.5.0

#environment variables
DBpig=00_data/00_databases/pig_gut/fna_files
DBown=00_data/04_bins_renamed
DBall=00_data/06_db_raw
DBdrep=00_data/06_db_drep
DBout=02_report/06_db_drep

#create directory
mkdir -p $DBall
mkdir -p $DBdrep
mkdir -p $DBout

# copy all bins to new directory

cp $DBpig/* $DBall
cp $DBown/* $DBall

# drep
dRep dereplicate $DBdrep -g $DBall/* -sa 0.95 -nc 0.3 -p 32 # tresholds as mgnify dbs

# move drep plots to reports folder

cp -r $DBdrep/figures $DBout
