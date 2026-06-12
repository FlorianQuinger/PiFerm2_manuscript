#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=16g
#SBATCH --job-name=047_rename
#SBATCH --output=logs/047_rename.out
#SBATCH --error=logs/047_rename.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics

module load devel/miniforge
conda activate bbmap_v39.17

#environment variables
BF=00_data/04_bins_filtered
BR=00_data/04_bins_renamed

#create directory
mkdir -p $BR

# rename bins
#i as counter
i=1
for FILE in $BF/*; do
	MAG=$(basename "$FILE")
 	MAG="${MAG/.fa/}"
 	MAGid=$(echo "$MAG" | grep -oE '^AP2_[0-9]{3}')
	MAGid=${MAGid/AP2_/}

 	newMAG="MGYGPF${MAGid}$i"

 	rename.sh in=$FILE out=$BR/$newMAG.fna prefix="$newMAG"

 	echo "renamed $MAG to $newMAG"

 	((i++))
done

# copy folder to report folder
cp -r $BR 02_report/
