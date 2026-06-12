#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=64g
#SBATCH --job-name=064_db_creation
#SBATCH --output=logs/064_db_creation.out
#SBATCH --error=logs/064_db_creation.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge

#environment variables
DB=02_report/06_db
DBmerge=$DB/pig_gut1.1_merge
DBdrep=$DB/pig_gut1.1_drep
DBown=$DB/pig_gut1.1_own
BP=00_data/06_bins_prokka
BEN=02_report/06_protein_eggnog
DREP=00_data/06_db_drep/dereplicated_genomes
#Odb=00_data/00_databases/pig_gut/original_db # replaced 
PGprokka=00_data/06_pig_gut_prokka
#Oegg=00_data/00_databases/pig_gut/eggNOG # replaced
PGegg=02_report/06_pig_gut_eggnog

#create directory
mkdir -p $DB
mkdir -p $DBmerge
mkdir -p $DBdrep
mkdir -p $DBown

# copy only dereplicated bins for drep DB

mkdir -p $DBdrep/original_db # for faa files
mkdir -p $DBdrep/eggNOG # for eggnog files

for i in $DREP/MGYG0*; do # copy bins from original db
	MAG=$(basename $i)
	MAG=${MAG/.fna/}
	cp $PGprokka/${MAG}.faa $DBdrep/original_db
	cp $PGegg/${MAG}.emapper.annotations $DBdrep/eggNOG/${MAG}_eggNOG.tsv
	echo "copied $i"
done

for i in $DREP/MGYGPF*; do # copy own bins
	MAG=$(basename $i)
	MAG=${MAG/.fna/}
	cp $BP/${MAG}.faa $DBdrep/original_db
	cp $BEN/${MAG}.emapper.annotations $DBdrep/eggNOG/${MAG}_eggNOG.tsv
	echo "copied $i"
done	

echo "finished creation of dereplicated database"