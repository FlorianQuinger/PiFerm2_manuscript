library(here)

source(here("0_general_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load pig_gut genomes-all_metadata.tsv file

metadata_pig_gut <- read_tsv("../../../Bioinformatics/MetaLab_MAG/databases/pig_gut/genomes-all_metadata.tsv")

# load reannotation with GTDBKtk 220

pig_gut_classification_ar <- read_tsv("../PiFerm2_metagenomics/02_report/06_pig_gut_gtdbtk/gtdbtk.ar53.summary.tsv") 
pig_gut_classification_ba <- read_tsv("../PiFerm2_metagenomics/02_report/06_pig_gut_gtdbtk/gtdbtk.bac120.summary.tsv") 

pig_gut_classification <- rbind(pig_gut_classification_ar, pig_gut_classification_ba) %>%
  dplyr::select(Genome = user_genome, classification)

# replace old taxonomy in metadata

unrelevant_genomes <- !(metadata_pig_gut$Genome %in% pig_gut_classification$Genome)

metadata_pig_gut <- metadata_pig_gut %>%
  inner_join(pig_gut_classification, by = "Genome") %>%
  mutate(Lineage = classification) %>%
  dplyr::select(-classification)

# create own metadata

# load bin stats

bin_stats <- read_tsv("../PiFerm2_metagenomics/02_report/05_assembly_stats/stats.txt") %>%
  mutate(bin = str_remove(filename, "/pfs/10/work/ho_yogau97-PiFerm2/PiFerm2_metagenomics/00_data/04_bins_renamed/"),
         bin = str_remove(bin, ".fna")) %>%
  dplyr::select(bin, length = contig_bp, contigs = n_contigs, N50 = ctg_L50, L50 = ctg_N50, gc = gc_avg) 

# load bin quality

# match renamed bins to previous bins with quality scores
bin_to_bin_renamed <- read_tsv("../PiFerm2_metagenomics/01_code/logs/047_rename.out", col_names = F) %>%
  mutate(clean = str_remove(X1, "renamed "), 
         clean = str_remove(clean, "to "),
         bin = str_extract(clean, "^[A-Za-z0-9_.]+"),
         bin_renamed = str_extract(clean, "[A-Z0-9]+$")) %>%
  dplyr::select(bin, bin_renamed)

bin_quality <- read_tsv("../PiFerm2_metagenomics/02_report/04_checkm2/quality_report.tsv") %>%
  inner_join(bin_to_bin_renamed, by = c("Name" = "bin")) %>%
  dplyr::select(bin = bin_renamed, Completeness, Contamination) 

# load bin classification

bin_classification_ar <- read_tsv("../PiFerm2_metagenomics/02_report/05_gtdbtk/gtdbtk.ar53.summary.tsv") 
bin_classification_ba <- read_tsv("../PiFerm2_metagenomics/02_report/05_gtdbtk/gtdbtk.bac120.summary.tsv") 

bin_classification <- rbind(bin_classification_ar, bin_classification_ba) %>%
  dplyr::select(bin = user_genome, classification)

# create own metadata

metadata_own <- bin_stats %>%
  inner_join(bin_quality, by = "bin") %>%
  inner_join(bin_classification, by = "bin") %>%
  mutate(Genome_type = "MAG", rRNA_5S = 0, rRNA_16S = 0, rRNA_23S = 0, tRNAs = 0, Genome_accession = NA, Species_rep = NA, Sample_accession = "PiFerm2", Study_accession = "PiFerm", Country = "Germany", Continent = "Europe", FTP_download = NA) %>%
  dplyr::select(Genome = bin, Genome_type, Length = length, N_contigs = contigs, N50, GC_content = gc, Completeness,
                Contamination, rRNA_5S, rRNA_16S, rRNA_23S, tRNAs, Genome_accession, Species_rep, 
                Lineage = classification, Sample_accession, Study_accession, Country, Continent, FTP_download)

# For database with dereplicated MAGs

drep <- str_remove(list.files("../../../Bioinformatics/MetaLab_MAG_1.1/db/pig_gut1.1_drep/original_db/"), ".faa$")

genomes_all_metadata_drep <- genomes_all_metadata_merge %>%
  filter(Genome %in% drep)

write_tsv(genomes_all_metadata_drep, "../../../Bioinformatics/MetaLab_MAG_1.1/db/pig_gut1.1_drep/genomes-all_metadata.tsv")

