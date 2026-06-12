#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=300g
#SBATCH --job-name=023_bbmap_db
#SBATCH --output=logs/023_bbmap_db.out
#SBATCH --error=logs/023_bbmap_db.err
cd /pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/00_databases

module load devel/miniforge
conda activate bbmap_v39.17

# use script to split chromosomes
for i in GCF_*; do
 OUT="split_${i}"
 zcat $i | awk -v max_length=450000000 -f ../../01_code/023_bbmap_db_split.awk > $OUT
done

# concatenate all output files

cat split_GCF* > reference_split.fna
gzip reference_split.fna

#converst SILVA to DNA
esl-reformat -d -o SILVA_dna.fasta fasta SILVA_138.2_Ref_NR99.fasta.gz
gzip SILVA_dna.fasta

# shred ribosomal sequences
shred.sh in=SILVA_dna.fasta.gz out=SILVA_shredded.fa.gz length=80 minlength=70 overlap=40

echo "shredding done"

# map to reference genomes
bbmap.sh ref=reference_split.fna.gz in=SILVA_shredded.fa.gz outm=reference_mapped.sam minid=0.9 maxindel=2
#bbmap.sh ref=GCF_000003025.6_Sscrofa11.1_genomic.fna.gz in=SILVA_shredded.fa.gz outm=reference_mapped.sam minid=0.9 maxindel=2

echo "mapping done"

# mask reference genomes
bbmask.sh in=reference_split.fna.gz out=reference_masked.fa.gz entropy=0.7 sam=reference_mapped.sam

echo "masking done"
