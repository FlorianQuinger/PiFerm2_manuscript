library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

meta <- readRDS("clean/meta1.RDS")

################################
# Reads
###############################

# gtdb v220 taxonomy

bac120 <- read_tsv("data/bac120_taxonomy_r220.tsv", col_names = F)
ar53 <- read_tsv("data/ar53_taxonomy_r220.tsv", col_names = F)

gtdb_v220 <- rbind(bac120, ar53) %>%
  dplyr::select(-X1) %>%
  mutate(X2 = str_remove_all(X2, "[a-z]__")) %>%
  separate(X2, into = c("K", "P", "C", "O", "F", "G", "S"), sep = ";") %>%
  distinct()

saveRDS(gtdb_v220, "clean/3_gtdb_v220_taxonomy.RDS")

# taxonomy of reads by kraken2/bracken

## k2report files

k2reports <- list.files("../PiFerm2_metagenomics/02_report/02_kraken2/", pattern = "k2report$")

k2report_files <- list()
for (i in 1:length(k2reports)) {
  filename <- k2reports[i]
  sampleid <- str_remove(filename, ".k2report")
  file <- read_tsv(paste0("../PiFerm2_metagenomics/02_report/02_kraken2/", filename), col_names = F)[,c(1,2, 4,6)] %>%
    dplyr::rename(rel_abd = X1, reads = X2, rank = X4, name = X6) %>% # rel_abds do not sum up to 100% on each level
    add_column(sampleid = sampleid)
  k2report_files[[i]] <- file
}

unclassified_kraken <- bind_rows(k2report_files) %>%
  dplyr::select(sampleid, name, rank, rel_abd) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  filter(rank == "U") %>%
  inner_join(meta, by = "sampleid") %>%
  ggplot(aes(x = matrix, y = rel_abd)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom(width = 0.3) +
  labs(title = "Unclassified reads by Kraken2")

print(unclassified_kraken)
save_big("04_kraken2_unclassified_reads")

k2report_reads <- bind_rows(k2report_files) %>%
  dplyr::select(sampleid, name, rank, reads) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  pivot_wider(names_from = "sampleid", values_from = "reads") %>%
  pivot_longer(-c("name", "rank"), names_to = "sampleid", values_to = "reads") %>%
  mutate(reads = ifelse(is.na(reads), 0, reads)) %>% # replace NA with 0
  group_by(rank, name) %>%
  filter(sum(reads) > 0) %>% # remove rows with only 0
  ungroup()

k2report_reads %>%
  filter(rank == "S") %>%
  summarize(reads = sum(reads)) # different ranks have different number of total reads, probably due to exclusion of reads that were not classified up to a specific rank

saveRDS(k2report_reads, file = "clean/3_k2report_reads_long_reads.RDS")

k2_report_reads_filtered <- filter_frequency_and_abundance_ktable(k2report_reads)

saveRDS(k2_report_reads_filtered, "clean/3_k2report_reads_long_reads_filtered.RDS")

k2report_rel_abd <- bind_rows(k2report_files) %>%
  dplyr::select(sampleid, name, rank, rel_abd) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  pivot_wider(names_from = "sampleid", values_from = "rel_abd") %>%
  pivot_longer(-c("name", "rank"), names_to = "sampleid", values_to = "rel_abd") %>%
  mutate(rel_abd = ifelse(is.na(rel_abd), 0, rel_abd)) %>% # replace NA with 0
  group_by(rank, name) %>%
  filter(sum(rel_abd) > 0) %>% # remove rows with only 0
  ungroup() %>%
  recalculate_rel_abd(rank_column = "rank")

saveRDS(k2report_rel_abd, file = "clean/3_k2report_reads_long_rel_abd.RDS")

# using reads to recalculate rel_abd and harmonise both versions after filtering

k2report_rel_abd <- k2report_reads %>%
  group_by(sampleid, rank) %>%
  mutate(rel_abd = reads / sum(reads) * 100) %>%
  ungroup() %>%
  dplyr::select(-reads)

k2_report_rel_abd_filtered <- filter_frequency_and_abundance_ktable(k2report_rel_abd)

saveRDS(k2_report_rel_abd_filtered, "clean/3_k2report_reads_long_rel_abd_filtered.RDS")

## breport files

breports <- list.files("../PiFerm2_metagenomics/02_report/02_kraken2/", pattern = "breport$")

breport_files <- list()
for (i in 1:length(breports)) {
  filename <- breports[i]
  sampleid <- str_remove(filename, ".breport")
  file <- read_tsv(paste0("../PiFerm2_metagenomics/02_report/02_kraken2/", filename), col_names = F)[,c(1,2,4,6)] %>%
    dplyr::rename(rel_abd = X1, reads = X2, rank = X4, name = X6) %>%
    add_column(sampleid = sampleid) # rel_abds sum up to 100 % at each level
  breport_files[[i]] <- file
}

breport_reads <- bind_rows(breport_files) %>%
  dplyr::select(sampleid, name, rank, reads) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  pivot_wider(names_from = "sampleid", values_from = "reads") %>%
  pivot_longer(-c("name", "rank"), names_to = "sampleid", values_to = "reads") %>%
  mutate(reads = ifelse(is.na(reads), 0, reads)) %>% # replace NA with 0
  group_by(rank, name) %>%
  filter(sum(reads) > 0) %>% # remove rows with only 0
  ungroup()

saveRDS(breport_reads, file = "clean/3_breport_reads_long_reads.RDS")

breport_reads_filtered <- filter_frequency_and_abundance_ktable(breport_reads)

saveRDS(breport_reads_filtered, "clean/3_breport_reads_long_reads_filtered.RDS")

breport_rel_abd <- bind_rows(breport_files) %>%
  dplyr::select(sampleid, name, rank, rel_abd) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  pivot_wider(names_from = "sampleid", values_from = "rel_abd") %>%
  pivot_longer(-c("name", "rank"), names_to = "sampleid", values_to = "rel_abd") %>%
  mutate(rel_abd = ifelse(is.na(rel_abd), 0, rel_abd)) %>% # replace NA with 0
  group_by(rank, name) %>%
  filter(sum(rel_abd) > 0) %>% # remove rows with only 0
  ungroup() %>%
  recalculate_rel_abd(rank_column = "rank")

test <- breport_rel_abd %>%
  group_by(rank, sampleid) %>%
  summarise(sum = sum(rel_abd))

saveRDS(breport_rel_abd, file = "clean/3_breport_reads_long_rel_abd.RDS")

# using reads to recalculate rel_abd and harmonise both versions after filtering

breport_rel_abd <- breport_reads %>%
  group_by(sampleid, rank) %>%
  mutate(rel_abd = reads / sum(reads) * 100) %>%
  ungroup() %>%
  dplyr::select(-reads)

breport_rel_abd_filtered <- filter_frequency_and_abundance_ktable(breport_rel_abd)

saveRDS(breport_rel_abd_filtered, "clean/3_breport_reads_long_rel_abd_filtered.RDS")

## bracken files

brackens <- list.files("../PiFerm2_metagenomics/02_report/02_kraken2/", pattern = "bracken$")

bracken_files <- list()
for (i in 1:length(brackens)) {
  filename <- brackens[i]
  sampleid <- str_remove(filename, ".bracken")
  file <- read_tsv(paste0("../PiFerm2_metagenomics/02_report/02_kraken2/", filename)) %>%
    add_column(sampleid = sampleid)
  bracken_files[[i]] <- file
}

bracken_reads <- bind_rows(bracken_files) %>% # only species level
  dplyr::select(sampleid, name, reads = new_est_reads) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  pivot_wider(names_from="sampleid", values_from="reads") %>%
  pivot_longer(-name, names_to="sampleid", values_to="reads") %>%
  mutate(reads = ifelse(is.na(reads), 0, reads)) %>%
  group_by(name) %>%
  filter(sum(reads) > 0) %>% # remove rows with only 0
  ungroup()

saveRDS(bracken_reads, file = "clean/3_bracken_reads_long_reads.RDS")

bracken_reads_filtered <- filter_frequency_and_abundance(bracken_reads)

saveRDS(bracken_reads_filtered, "clean/3_bracken_reads_long_reads_filtered.RDS")

bracken_rel_abd <- bind_rows(bracken_files) %>% # only species level
  dplyr::select(sampleid, name, rel_abd = fraction_total_reads) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  pivot_wider(names_from="sampleid", values_from="rel_abd") %>%
  pivot_longer(-name, names_to="sampleid", values_to="rel_abd") %>%
  mutate(rel_abd = ifelse(is.na(rel_abd), 0, rel_abd)) %>%
  group_by(name) %>%
  filter(sum(rel_abd) > 0) %>% # remove rows with only 0
  ungroup() %>%
  recalculate_rel_abd()

saveRDS(bracken_rel_abd, file = "clean/3_bracken_reads_long_rel_abd.RDS")

# using reads to recalculate rel_abd and harmonise both versions after filtering

bracken_rel_abd <- bracken_reads %>%
  group_by(sampleid) %>%
  mutate(rel_abd = reads / sum(reads) * 100) %>%
  ungroup() %>%
  dplyr::select(-reads)

bracken_rel_abd_filtered <- filter_frequency_and_abundance(bracken_rel_abd)

saveRDS(bracken_rel_abd_filtered, "clean/3_bracken_reads_long_rel_abd_filtered.RDS")

#############################################
# Contigs
#############################################

# function of contigs

## quantifications of genes

coverm_genes <- list.files("../PiFerm2_metagenomics/02_report/03_contigs_coverm_eggnog/")

coverm_genes_files <- list()
for (i in 1:length(coverm_genes)) {
  filename <- coverm_genes[i]
  sampleid <- str_remove(filename, ".txt")
  file <- read_tsv(paste0("../PiFerm2_metagenomics/02_report/03_contigs_coverm_eggnog/", filename)) %>%
    add_column(sampleid = sampleid)
  colnames(file) <- c("gene", "rel_abd", "reads", "sampleid")
  coverm_genes_files[[i]] <- file
}

coverm_gene <- bind_rows(coverm_genes_files) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_"))

unmapped_genes <- coverm_gene %>%
  filter(gene == "unmapped") %>%
  inner_join(meta, by = "sampleid") %>%
  ggplot(aes(x = matrix, y = rel_abd)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom(width = 0.3) +
  labs(title = "Unmapped reads to contig genes")

print(unmapped_genes)
save_big("04_coverm_unmapped_genes")

## eggNOG

eggnogs <- list.files("../PiFerm2_metagenomics/02_report/03_eggnog/", pattern = ".emapper.annotations")

eggnog_files <- list()
for (i in 1:length(eggnogs)) {
  filename <- eggnogs[i]
  sampleid <- str_remove(filename, ".emapper.annotations")
  file <- read_tsv(paste0("../PiFerm2_metagenomics/02_report/03_eggnog/", filename)) %>%
    add_column(sampleid = sampleid)
  eggnog_files[[i]] <- file
}

eggnog_meta <- bind_rows(eggnog_files) %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  dplyr::select(-c("#query", sampleid, evalue, score, sampleid)) %>% 
  filter(str_detect(eggnog_ogs, "Bacteria|Archaea")) %>% # only select bacterial and archaeal functions (all of them)
  distinct() %>%
  mutate(cog_accession = str_extract_all(eggnog_ogs, "COG\\d+@1"), # extracting root COGs
         cog_accession = str_remove_all(cog_accession, "@1"), # remove ndicators for root anc coerce
         cog_accession = ifelse(cog_accession == "character(0)", NA, str_remove_all(cog_accession, '[c()"\\s]'))) # insert NA if no COG exists and remove vector characters

eggnog <- bind_rows(eggnog_files) %>% 
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_")) %>%
  dplyr::rename(gene = "#query") %>%
  filter(str_detect(eggnog_ogs, "Bacteria|Archaea")) %>% # only select bacterial and archaeal functions (all of them)
  left_join(coverm_gene, by = c("sampleid", "gene")) %>%
  mutate(cog_accession = str_extract_all(eggnog_ogs, "COG\\d+@1"), # extracting root COGs
         cog_accession = str_remove_all(cog_accession, "@1"), # remove ndicators for root anc coerce
         cog_accession = ifelse(cog_accession == "character(0)", NA, str_remove_all(cog_accession, '[c()"\\s]'))) # insert NA if no COG exists and remove vector characters
  
saveRDS(eggnog, file = "clean/3_eggnog_contigs_long_raw.RDS")

# plot classified genes

eggnog_classified <- eggnog %>%
  group_by(sampleid) %>%
  summarise(rel_abd = sum(rel_abd)) %>%
  inner_join(meta, by = "sampleid")

ggplot(eggnog_classified, aes(x = matrix, y = rel_abd)) +
  geom_boxplot() +
  geom_quasirandom(aes(color = animal, shape = diet), size = 2) +
  scale_color_manual(values = colors)

# combine gene function with taxonomic information about contigs

eggnog_mmseqs2 <- eggnog %>%
  mutate(contig = str_extract(gene, "k\\d+_\\d+")) %>%
  left_join(mmseqs2_meta, by = c("sampleid", "contig"))

saveRDS(eggnog_mmseqs2, "clean/3_eggnog_mmseqs2_contigs_combined.RDS")

# merge 

# rel_abd and reads

eggnog_reads <- eggnog %>%
  dplyr::select(sampleid, seed_ortholog, reads) %>%
  group_by(sampleid, seed_ortholog) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  pivot_wider(names_from = sampleid, values_from = reads) %>%
  pivot_longer(-seed_ortholog, names_to = "sampleid", values_to = "reads") %>%
  mutate(reads = ifelse(is.na(reads), 0, reads)) %>%
  inner_join(eggnog_meta, by = "seed_ortholog")

saveRDS(eggnog_reads, file = "clean/3_eggnog_contigs_long_reads.RDS")

eggnog_reads_filtered <- filter_frequency_and_abundance(eggnog_reads, frequency_cutoff = 0.10, 
                                                        abundance_cutoff = 1e-4) # lower thresholds for functions

saveRDS(eggnog_reads_filtered, file = "clean/3_eggnog_contigs_long_reads_filtered.RDS")

eggnog_rel_abd <- eggnog %>%
  dplyr::select(sampleid, seed_ortholog, rel_abd) %>%
  group_by(sampleid, seed_ortholog) %>%
  summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = sampleid, values_from = rel_abd) %>%
  pivot_longer(-seed_ortholog, names_to = "sampleid", values_to = "rel_abd") %>%
  mutate(rel_abd = ifelse(is.na(rel_abd), 0, rel_abd)) %>%
  inner_join(eggnog_meta, by = "seed_ortholog") %>%
  recalculate_rel_abd()

saveRDS(eggnog_rel_abd, file = "clean/3_eggnog_contigs_long_rel_abd.RDS")

eggnog_rel_abd_filtered <- filter_frequency_and_abundance(eggnog_rel_abd, frequency_cutoff = 0.10, 
                                                          abundance_cutoff = 1e-4) 

saveRDS(eggnog_rel_abd_filtered, file = "clean/3_eggnog_contigs_long_rel_abd_filtered.RDS")

# to KEGG

kegg_reads <- calculate_kegg_abundance(eggnog_reads, abundance_column = "reads")
saveRDS(kegg_reads, "clean/3_kegg_contigs_long_reads.RDS")

#kegg_reads_filtered <- filter_frequency_and_abundance(kegg_reads) # alternative method for filtered frame
kegg_reads_filtered <- calculate_kegg_abundance(eggnog_reads_filtered, abundance_column = "reads")
saveRDS(kegg_reads_filtered, "clean/3_kegg_contigs_long_reads_filtered.RDS")

kegg_rel_abd <- calculate_kegg_abundance(eggnog_rel_abd, abundance_column = "rel_abd") %>%
  recalculate_rel_abd()
saveRDS(kegg_rel_abd, "clean/3_kegg_contigs_long_rel_abd.RDS") #TODO

#kegg_rel_abd_filtered <- filter_frequency_and_abundance(kegg_rel_abd)
kegg_rel_abd_filtered <- calculate_kegg_abundance(eggnog_rel_abd_filtered, abundance_column = "rel_abd") %>%
  recalculate_rel_abd()
saveRDS(kegg_rel_abd_filtered, "clean/3_kegg_contigs_long_rel_abd_filtered.RDS")

# to Go

go_reads <- calculate_go_abundance(eggnog_reads, abundance_column = "reads")
saveRDS(go_reads, "clean/3_go_contigs_long_reads.RDS")

#go_reads_filtered <- filter_frequency_and_abundance(go_reads)
go_reads_filtered <- calculate_go_abundance(eggnog_reads_filtered, abundance_column = "reads")
saveRDS(go_reads_filtered, "clean/3_go_contigs_long_reads_filtered.RDS")

go_rel_abd <- calculate_go_abundance(eggnog_rel_abd, abundance_column = "rel_abd") %>%
  recalculate_rel_abd()
saveRDS(go_rel_abd, "clean/3_go_contigs_long_rel_abd.RDS")

#go_rel_abd_filtered <- filter_frequency_and_abundance(go_rel_abd)
go_rel_abd_filtered <- calculate_go_abundance(eggnog_rel_abd_filtered, abundance_column = "rel_abd")
saveRDS(go_rel_abd_filtered, "clean/3_go_contigs_long_rel_abd_filtered.RDS")

# to COG

cog_reads <- calculate_cog_abundance(eggnog_reads, abundance_column = "reads")
saveRDS(cog_reads, "clean/3_cog_contigs_long_reads.RDS")

#cog_reads_filtered <- filter_frequency_and_abundance(cog_reads)
cog_reads_filtered <- calculate_cog_abundance(eggnog_reads_filtered, abundance_column = "reads")
saveRDS(cog_reads_filtered, "clean/3_cog_contigs_long_reads_filtered.RDS")

cog_rel_abd <- calculate_cog_abundance(eggnog_rel_abd, abundance_column = "rel_abd") %>%
  recalculate_rel_abd()
saveRDS(cog_rel_abd, "clean/3_cog_contigs_long_rel_abd.RDS")

#cog_rel_abd_filtered <- filter_frequency_and_abundance(cog_rel_abd)
cog_rel_abd_filtered <- calculate_cog_abundance(eggnog_rel_abd_filtered, abundance_column = "rel_abd")
saveRDS(cog_rel_abd_filtered, "clean/3_cog_contigs_long_rel_abd_filtered.RDS")

# to CAZY

cazy_reads <- calculate_cazy_abundance(eggnog_reads, abundance_column = "reads")
saveRDS(cazy_reads, "clean/3_cazy_contigs_long_reads.RDS")

#cazy_reads_filtered <- filter_frequency_and_abundance(cazy_reads)
cazy_reads_filtered <- calculate_cazy_abundance(eggnog_reads_filtered, abundance_column = "reads")
saveRDS(cazy_reads_filtered, "clean/3_cazy_contigs_long_reads_filtered.RDS")

cazy_rel_abd <- calculate_cazy_abundance(eggnog_rel_abd, abundance_column = "rel_abd") %>%
  recalculate_rel_abd()
saveRDS(cazy_rel_abd, "clean/3_cazy_contigs_long_rel_abd.RDS")

#cazy_rel_abd_filtered <- filter_frequency_and_abundance(cazy_rel_abd)
cazy_rel_abd_filtered <- calculate_cazy_abundance(eggnog_rel_abd_filtered, abundance_column = "rel_abd")
saveRDS(cazy_rel_abd_filtered, "clean/3_cazy_contigs_long_rel_abd_filtered.RDS")

