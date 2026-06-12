library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta

meta <- readRDS("clean/meta1.RDS") %>% filter(sampleid != "103")

# load data

#protein_func_tax <- readRDS("clean/4_function_taxa_combined_long_norm_imp_rel_abd_filtered.RDS")
protein_func_tax <- readRDS("clean/4_function_taxa_combined_long_raw_rel_abd_filtered.RDS")

# load taxonomy files

#taxonomy <- readRDS("clean/4_taxonomy_long_norm_imp_rel_abd_filtered.RDS") 
taxonomy <- readRDS("clean/4_taxonomy_long_raw_rel_abd_filtered.RDS")

# create protein-kegg/cog-genus table

pro_gen_cog <- protein_func_tax %>%
  mutate(kegg_cog = ifelse(kegg_ko == "-", cog_accession, kegg_ko)) %>%
  dplyr::select(proteinid, G, kegg_cog, sampleid, rel_abd) %>%
  mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
  separate_rows(kegg_cog, sep = ",") %>%
  group_by(id) %>%
  mutate(factor = 1 / length(kegg_cog)) %>% # calculating factor
  ungroup() %>%
  mutate(rel_abd = rel_abd * factor) %>%
  dplyr::select(-id, -factor) %>%
  mutate(kegg_cog = str_remove(kegg_cog, "ko:")) 

genus_table <- taxonomy %>%
  filter(rank == "G") %>%
  dplyr::select(-rank)



nfr <- calculate_functional_redundancy(pro_gen_cog = pro_gen_cog, genus_table = genus_table, meta = meta) %>%
  as_tibble(rownames = "sampleid") %>%
  inner_join(meta, by = "sampleid")

ggplot(nfr, aes(x = diet, y = nFR)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom(aes(shape = animal)) +
  facet_wrap(vars(matrix), scales = "free_y")

#comparison

nfr_il <- select_matrix(nfr, "ileal digesta")
nfr_fa <- select_matrix(nfr, "faeces")

comp_nfr_il <- combined_comparison(select_response(nfr_il, "nFR"), transformation = "test") 
comp_nfr_fa <- combined_comparison(select_response(nfr_fa, "nFR"), transformation = "test")

table_nfr_il <- create_results_table(select_response(nfr_il, "nFR"), comp_nfr_il, response = "nFR_il", digits = 2)
table_nfr_fa <- create_results_table(select_response(nfr_fa, "nFR"), comp_nfr_fa, response = "nFR_fa", digits = 2)

table_nfr <- table_nfr_il %>%
  inner_join(table_nfr_fa, by = "diet")

write_tsv(table_nfr, "tables/49_nfr.txt")