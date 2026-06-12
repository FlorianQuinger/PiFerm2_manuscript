#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=500g
#SBATCH --job-name=066_pig_gut_gtdbtk
#SBATCH --output=logs/066_pig_gut_gtdbtk.out
#SBATCH --error=logs/066_pig_gut_gtdbtk.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate gtdbtk_v2.4.0

#environment variable for database
conda env config vars set GTDBTK_DATA_PATH="00_data/00_databases/gtdbtk_r220"

#environment variables
DB=00_data/00_databases/gtdbtk_r220
PGF=00_data/00_databases/pig_gut/fna_files
TG=00_data/06_pig_gut_taxonomy_gtdbtk
OUT=02_report/06_pig_gut_gtdbtk

#create directory
mkdir -p $DB
mkdir -p $TG
mkdir -p $TG/iqtree
mkdir -p $OUT

# download database
#wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
#tar -xvzf gtdbtk_data.tar.gz -C $DB --strip 1 > /dev/null
#rm gtdbtk_data.tar.gz
conda env config vars set GTDBTK_DATA_PATH=$DB

# reactivate environment
conda activate gtdbtk_v2.4.0

#run gtdb classify workflow on all filtered MAGs
gtdbtk classify_wf --genome_dir $PGF --out_dir $TG --cpus 32 -x fna --skip_ani_screen #--mash_db $DB/gtdbtk_data.msh 

# copy output files to reports

cp $TG/*.summary.tsv $OUT/

# creating tree using iqtree and protein alignments
conda activate iqtree_v2.4.0

iqtree -s $TG/align/gtdbtk.bac120.user_msa.fasta.gz -nt 32 -pre $TG/iqtree/bac120
iqtree -s $TG/align/gtdbtk.ar53.user_msa.fasta.gz -nt 32 -pre $TG/iqtree/ar53

# copy tree files to output folder

cp $TG/iqtree/*.treefile $OUT/
