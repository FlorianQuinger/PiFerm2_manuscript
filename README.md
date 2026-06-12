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
- **064_db_creation.sh** - Copies predicted proteins and functional annotation files to create database

### statistics
The names of the diets may be different compared to the manuscript:
- diet1 = spring field pea = SP
- diet2 = winter field pea 1 = WP1
- diet3 = winter field pea 2 = WP2
- diet4 = spring forage pea = SFP

#### General files
- **00_experimental_design.R** - Randomization of the double Latin square design
- **0_general_functions.R** - Functions that may be reused across the whole analysis
- **30_omics_functions.R** - Functions used in script for metagenomics and metaproteomics
 
#### Data cleaning
- **00_experimental_design.R** - Randomization of the double Latin square design
- **01_data_cleaning.R** - Reads in metabolomics as tsv files, cleans column names, calculates concentration in DM, and save data as RDS file
- **02_data_cleaning_metaproteomics.R** - Reading in metaproteomics output files from MetaLab MAG, cleaning, normalization, imputation, and generation of different output files for taxonomic and functional analysis
- **03_database_creation_metagenomes.R** - Helper script to create taxonomy file compatible with MetaLab MAG for dereplicated MAG database
- **04_data_cleaning_metagenomics.R** - Reading in output files from metagenomics pipeline, cleaning, and generation of different output files

#### Enzyme activity
- **10_enzyme_functions.R** - Functions used in scripts for enzyme analysis
- **16_enzyme_dry_matter.R** - Calculating enzyme activity per DM
- **17_enzyme_analysis.R** - Statistical analysis of enzyme activities

#### Nutrient digestibility
- **20_nutrition_functions.R** - Functions usd in scripts for analysis of nutrient digestibility
- **21_nutrition_analysis.R** - Statistical analysis of data retrieved from nutrient analyses

#### Metagenomics
- **32_metagenomics_reads_taxonomy.R** - Visualization and statistical analysis of taxonomy based on reads
- **33_metagenomics_contigs_tax_func.R** - Visualization and statistical analysis of functional data based on contigs

#### Metaproteomics
- **44_metaproteomics_basic.R** - Visualization of metrics and micro-host ratio for metaproteomics
- **45_metaproteomics_ordination.R** - Beta diversity of metaproteomics data
- **46_metaproteomics_dea.R** - Differential abundance analysis of functions and proteins using edgeR
- **47_metaproteomics_daa.R** - Alpha diversity and differential abundance analysis of taxa using ANCOM-BC2
- **49_metaproteomics_nfr.R** - Calculation and statistical analysis of functional redundancy

#### Metabolomics
- **52_metabolomics_reml.R** - Statistical analysis of metabolite concentrations

#### Multiomics
- **70_multiomics_functions.R** - Functions that are used for multiomics analysis
- **70_multiomics_data_prep.R** - Preprocessing of the input data for DIABLO
- **73_multiomics_diablo_ileum.R** - DIABLO models for ileal digesta
- **74_multiomics_diablo_faeces.R** - DIABLO models for feces

#### Other
- **92_figures_JAS.R** - Creation of figures for publication
