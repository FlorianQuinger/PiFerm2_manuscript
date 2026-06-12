#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=046_checkm2
#SBATCH --output=logs/046_checkm2.out
#SBATCH --error=logs/046_checkm2.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate checkm2_v1.0.2

#environment variables
DB=00_data/00_databases
BDR=00_data/04_bins_drep/dereplicated_genomes
CM=02_report/04_checkm2
BF=00_data/04_bins_filtered

#create directory
mkdir -p $CM
mkdir -p $BF

# download checkm2 database
checkm2 database --download --path $DB

#checkm2 for all fasta in drep output
checkm2 predict --threads 32 --input $BDR/*.fa --output-directory $CM --database_path $DB/CheckM2_database/uniref100.KO.1.dmnd

# Filter bins for 50% completeness and 5% contamination
awk -F'\t' -v min_comp=50 -v max_cont=5 'NR>1 {
	file = $1
	completeness = $2
	contamination = $3

	if (completeness >= min_comp && contamination <= max_cont) {
		print file;
	}
}' $CM/quality_report.tsv | while read -r BIN; do
	# move files
	if [ -f "$BDR/$BIN".fa ]; then
		cp "$BDR/$BIN".fa "$BF/"
		echo "copied $BDR/$BIN"
	else
		echo "No $BIN found"
	fi
done

# copy folder to report folder
cp -r $BF 02_report/
