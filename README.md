# PiFerm2_manuscript

## Overview
This repository contains the code for the paper X. 
The content is separated in two directories. The directory "metagenomics" contains the code used for processing the raw shotgun sequencing files and the generation of the MAG database for metaproteomics. The code was ran on the bwForCluster BinAC2. The directory "statistics" contains all the analysis performed in R, including data wrangling, statistical analysis, and visualization.

## Content

### metagenomics

#### Data Validation & Preparation
- **00_checksum.sh** - Validates integrity of uploaded tar files using MD5 checksums
- **011_unzip.sh** - Extracts tar archives containing raw sequencing data
- **012_checksum.sh** - Verifies checksums of extracted sequencing files
- **013_move_reads.sh** - Consolidates sequencing files from subdirectories
- **014_concat_reads.sh** - Concatenates forward and backward sequencing files of paired-end reads into single files

### Quality Control & Read Processing
- **021_fastqc.sh** - Performs initial quality assessment using FastQC and MultiQC
- **022_bbduk.sh** - Trims adapters, filters contaminants (phix/artifacts), and performs quality trimming
- **023_bbmap_db.sh** - Prepares masked reference database for host genome removal
- **023_bbmap_db_split.awk** - Utility script to split large FASTA sequences into manageable chunks
- **023_bbmap.sh** - Maps reads to masked pig genome to remove host-derived sequences

### Taxonomic Annotation
- **024_kraken2_unzip.sh** - Downloads and extracts Kraken2 GTDB database
- **024_kraken2.sh** - Classifies reads using Kraken2 and estimates abundance with Bracken

### Assembly & Functional Annotation
- **031_megahit.sh** - Performs metagenomic assembly with MEGAHIT
- **0312_filter_contigs.sh** - Filters assembled contigs (minimum 2500 bp length)
- **032_metaquast.sh** - Evaluates assembly quality metrics 
- **037_eggnog.sh** - Functional annotation of contigs using eggNOG mapper
- **039_coverm_eggnog.sh** - Maps reads to predicted genes from contigs to calculate abundance

### Metagenomic Binning
- **033_coverm.sh** - Calculates coverage for metabat2
- **034_bowtie2.sh** - Generates filtered and sorted BAM alignment files for SemiBin2
- **041_metabat2.sh** - Bins contigs using MetaBAT2 
- **042_semibin2.sh** - Performs binning with SemiBin2
- **043_comebin_cuda.sh** - Binning using COMEBin
- **044_das_tool.sh** - Combines the binning results from the three binners and selects the best bins.
- **045_drep.sh** - Dereplicates bins (98% ANI threshold)
- **046_checkm2.sh** - Assesses genome quality and filters bins (≥50% completeness, ≤5% contamination)
- **047_rename_bins.sh** - Renames filtered bins
- **054_assembly_stats.sh** - Computes comprehensive assembly statistics for final bins

### Database Creation for Metaproteomics
- **051_gtdbtk.sh** - Taxonomically classifies own bins 
- **061_database_drep.sh** - Creates dereplicated database from MGnify pig-gut 1.0 database and own bins
- **062_prokka.sh** - Performs gene prediction on own bins
- **063_protein_eggnog.sh** - Functionally annotates predicted proteins using eggNOG mapper
- **065_pig_gut_prokka.sh** - Gene prediction on published MGnify pig-gut 1.0 bins
- **066_pig_gut_gtdbtk.sh** - Taxonomic classification of MGnify pig-gut 1.0 bins
- **067_pig_gut_eggnog.sh** - Functional annotation of MGnify pig-gut 1.0 bins
- **064_db_creation.sh** - Copies predicted proteins and functional annotation files to create different databases
