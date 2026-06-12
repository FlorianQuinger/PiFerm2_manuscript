#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --time=144:00:00
#SBATCH --mem=200g
#SBATCH --job-name=063_eggnog
#SBATCH --output=logs/063_eggnog.out
#SBATCH --error=logs/063_eggnog.err

module load devel/miniforge
conda activate eggnog_v2.1.12

cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/ 

#environment variables
BP=00_data/06_bins_prokka
OUT=02_report/06_protein_eggnog

#create directory
mkdir -p $OUT

# add to path
export PATH=~/.conda/envs/eggnog_v2.1.12:~/.conda/envs/eggnog_v2.1.12/bin:"$PATH"

# download eggnog databases
export EGGNOG_DATA_DIR=$TMPDIR
download_eggnog_data.py -y

# run emapper per sample
for FILE in $BP/*.faa; do 
 # sampleid
 MAG=$(basename "$FILE")
 MAG="${MAG/.faa/}" 
 
 if [ -f "$OUT/$MAG.emapper.hits" ]; then #check if output files already generated
	 echo "already processed $MAG"
 else
 	 echo "start emapper for $MAG"

 	 emapper.py -i $FILE --itype proteins -m diamond --no_file_comments --cpu 128 --dbmem --tax_scope "prokaryota_broad" --temp_dir $TMPDIR --output_dir $OUT -o $MAG # as in embl ebi pipeline
 	 echo "processed $MAG"
 fi
done
