#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=48:00:00
#SBATCH --mem=200g
#SBATCH --job-name=037_eggnog
#SBATCH --output=logs/037_eggnog.out
#SBATCH --error=logs/037_eggnog.err

module load devel/miniforge
conda activate eggnog_v2.1.12

cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/ 

#environment variables
FC=00_data/03_contigs_filtered
OUT=02_report/03_eggnog

#create directory
mkdir -p $OUT

# add to path
export PATH=~/.conda/envs/eggnog_v2.1.12:~/.conda/envs/eggnog_v2.1.12/bin:"$PATH"

# download eggnog databases
export EGGNOG_DATA_DIR=$TMPDIR
#download_eggnog_data.py -y

# run emapper per sample
for FOLDER in $FC/*; do 
 # sampleid
 sampleid=$(basename "$FOLDER")
 
 echo "start emapper for $FOLDER"

 emapper.py -i $FOLDER/final.contigs.fa --itype genome -m diamond --no_file_comments --cpu 0 --dbmem --tax_scope "prokaryota_broad" --temp_dir $TMPDIR --output_dir $OUT -o $sampleid # as in embl ebi pipeline
 #emapper.py -m diamond --itype genome -i $FOLDER/final.contigs.fa -o test

 echo "processed $sampleid"
done
