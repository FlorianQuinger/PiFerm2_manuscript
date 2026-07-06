library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta

meta <- readRDS("clean/meta1.RDS") %>% filter(sampleid != "103")

# Metagenomics

breport_rel_abd_filtered <- readRDS("clean/3_breport_reads_long_rel_abd_filtered.RDS")%>% filter(sampleid != "103")
breport_reads_filtered <- readRDS("clean/3_breport_reads_long_reads_filtered.RDS")%>% filter(sampleid != "103")


kegg_rel_abd_filtered <- readRDS("clean/3_kegg_contigs_long_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
kegg_reads_filtered <- readRDS("clean/3_kegg_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")


# create taxa tables for ileum

table_g_g_il <- filter_ktable(ktable = breport_rel_abd_filtered, meta = meta, selected_rank = "G", 
                            selected_matrix = "ileal digesta") %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(Genus = name, SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_g_g_il, "tables/93_table_g_g_il.txt")

table_g_s_il <- filter_ktable(ktable = breport_rel_abd_filtered, meta = meta, selected_rank = "S", 
                            selected_matrix = "ileal digesta") %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(Species = name, SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_g_s_il, "tables/93_table_g_s_il.txt")

# create taxa tables for faeces

table_g_g_fa <- filter_ktable(ktable = breport_rel_abd_filtered, meta = meta, selected_rank = "G", 
                            selected_matrix = "faeces") %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(Genus = name, SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_g_g_fa, "tables/93_table_g_g_fa.txt")

table_g_s_fa <- filter_ktable(ktable = breport_rel_abd_filtered, meta = meta, selected_rank = "S", 
                            selected_matrix = "faeces") %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(Species = name, SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_g_s_fa, "tables/93_table_g_s_fa.txt")

# create KEGG tables for ileum

table_g_kegg_il <- kegg_rel_abd_filtered %>%
  filter_ileum() %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_g_kegg_il, "tables/93_table_g_kegg_il.txt")

# create KEGG tables for faeces

table_g_kegg_fa <- kegg_rel_abd_filtered %>%
  filter_faeces() %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_g_kegg_fa, "tables/93_table_g_kegg_fa.txt")

# results of taxa differential abundance analysis

il_ancom_G <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "G", 
                                       selected_matrix = "ileal digesta")

table_g_il_ancom_g <- do.call(rbind, il_ancom_G$res_pair_list)
write_tsv(table_g_il_ancom_g, "tables/93_table_g_il_ancom_g.txt")

il_ancom_S <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "S", 
                                       selected_matrix = "ileal digesta")

table_g_il_ancom_s <- do.call(rbind, il_ancom_S$res_pair_list)
write_tsv(table_g_il_ancom_s, "tables/93_table_g_il_ancom_s.txt")

fa_ancom_G <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "G", 
                                       selected_matrix = "faeces")

table_g_fa_ancom_g <- do.call(rbind, fa_ancom_G$res_pair_list)
write_tsv(table_g_fa_ancom_g, "tables/93_table_g_fa_ancom_g.txt")

fa_ancom_S <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "S", 
                                       selected_matrix = "faeces")

table_g_fa_ancom_s <- do.call(rbind, fa_ancom_S$res_pair_list)
write_tsv(table_g_fa_ancom_s, "tables/93_table_g_fa_ancom_s.txt")

# results of kegg differential abundance analysis

# differential abundance of keggs

il_kegg <- perform_edger_and_plot(kegg_reads_filtered, selected_matrix = "ileum", save = save,
                                  abundance_column = "reads", save_name = "kegg_contigs")

table_g_il_edger_kegg <- rbind(il_kegg$`diet2-diet1`$table_sig, il_kegg$`diet3-diet1`$table_sig, 
                               il_kegg$`diet4-diet1`$table_sig, il_kegg$`diet3-diet2`$table_sig, 
                               il_kegg$`diet4-diet2`$table_sig, il_kegg$`diet4-diet3`$table_sig)
write_tsv(table_g_il_edger_kegg, "tables/93_table_g_il_edger_kegg.txt")

fa_kegg <- perform_edger_and_plot(kegg_reads_filtered, selected_matrix = "faeces", save = save,
                                  abundance_column = "reads", save_name = "kegg_contigs")

table_g_fa_edger_kegg <- rbind(fa_kegg$`diet2-diet1`$table_sig, fa_kegg$`diet3-diet1`$table_sig, 
                               fa_kegg$`diet4-diet1`$table_sig, fa_kegg$`diet3-diet2`$table_sig, 
                               fa_kegg$`diet4-diet2`$table_sig, fa_kegg$`diet4-diet3`$table_sig)
write_tsv(table_g_fa_edger_kegg, "tables/93_table_g_fa_edger_kegg.txt")


# Metaproteomics

taxonomy_norm_imp_rel_abd <- readRDS("clean/4_taxonomy_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
taxonomy_norm_imp_intensity <- readRDS("clean/4_taxonomy_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")

kegg_norm_imp <- readRDS("clean/4_kegg_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
kegg_norm_imp_rel_abd <- readRDS("clean/4_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

host_kegg_norm_imp <- readRDS("clean/4_host_kegg_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
host_kegg_norm_imp_rel_abd <- readRDS("clean/4_host_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

# create lists for taxa in ileal digesta

table_p_g_il <- filter_ktable(ktable = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", 
                              selected_matrix = "ileal digesta") %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(Genus = name, SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_p_g_il, "tables/93_table_p_g_il.txt")

# create lists for taxa in faeces

table_p_g_fa <- filter_ktable(ktable = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", 
                              selected_matrix = "faeces") %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(Genus = name, SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_p_g_fa, "tables/93_table_p_g_fa.txt")

# tables for kegg abundances

table_p_kegg_il <- kegg_norm_imp_rel_abd %>%
  filter_ileum() %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_p_kegg_il, "tables/93_table_p_kegg_il.txt")

table_p_kegg_fa <- kegg_norm_imp_rel_abd %>%
  filter_faeces() %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_p_kegg_fa, "tables/93_table_p_kegg_fa.txt")

# tables for host kegg abundances

table_p_host_kegg_il <- host_kegg_norm_imp_rel_abd %>%
  filter_ileum() %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_p_host_kegg_il, "tables/93_table_p_host_kegg_il.txt")

table_p_host_kegg_fa <- host_kegg_norm_imp_rel_abd %>%
  filter_faeces() %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  pivot_wider(names_from = "diet", values_from = "rel_abd") %>%
  rename(SP = `1`, WP1 = `2`, WP2 = `3`, SFP = `4`)
write_tsv(table_p_host_kegg_fa, "tables/93_table_p_host_kegg_fa.txt")

# differential abundance of taxa

il_ancom_G <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "G", 
                                       selected_matrix = "ileal digesta")

table_p_il_ancom_g <- do.call(rbind, il_ancom_G$res_pair_list)
write_tsv(table_p_il_ancom_g, "tables/93_table_p_il_ancom_g.txt")

fa_ancom_G <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "G", 
                                       selected_matrix = "faeces")

table_p_fa_ancom_g <- do.call(rbind, fa_ancom_G$res_pair_list)
write_tsv(table_p_fa_ancom_g, "tables/93_table_p_fa_ancom_g.txt")

# differential abundance of keggs

kegg_il_results <- perform_edger_and_plot(kegg_norm_imp, selected_matrix = "ileum", save = F,
                                          abundance_column = "intensity", save_name = "kegg_micro")

table_p_il_edger_kegg <- rbind(kegg_il_results$`diet2-diet1`$table_sig, kegg_il_results$`diet3-diet1`$table_sig, 
                               kegg_il_results$`diet4-diet1`$table_sig, kegg_il_results$`diet3-diet2`$table_sig, 
                               kegg_il_results$`diet4-diet2`$table_sig, kegg_il_results$`diet4-diet3`$table_sig)
write_tsv(table_p_il_edger_kegg, "tables/93_table_p_il_edger_kegg.txt")

kegg_fa_results <- perform_edger_and_plot(kegg_norm_imp, selected_matrix = "faeces", save = F,
                                          abundance_column = "intensity", save_name = "kegg_micro")

table_p_fa_edger_kegg <- rbind(kegg_fa_results$`diet2-diet1`$table_sig, kegg_fa_results$`diet3-diet1`$table_sig, 
                               kegg_fa_results$`diet4-diet1`$table_sig, kegg_fa_results$`diet3-diet2`$table_sig, 
                               kegg_fa_results$`diet4-diet2`$table_sig, kegg_fa_results$`diet4-diet3`$table_sig)
write_tsv(table_p_fa_edger_kegg, "tables/93_table_p_fa_edger_kegg.txt")

# differential abundance of host keggs

kegg_host_il_results <- perform_edger_and_plot(host_kegg_norm_imp, selected_matrix = "ileum", save = F,
                                               abundance_column = "intensity", save_name = "kegg_host")

table_p_host_il_edger_kegg <- rbind(kegg_host_il_results$`diet2-diet1`$table_sig, kegg_host_il_results$`diet3-diet1`$table_sig, 
                               kegg_host_il_results$`diet4-diet1`$table_sig, kegg_host_il_results$`diet3-diet2`$table_sig, 
                               kegg_host_il_results$`diet4-diet2`$table_sig, kegg_host_il_results$`diet4-diet3`$table_sig)
write_tsv(table_p_host_il_edger_kegg, "tables/93_table_p_host_il_edger_kegg.txt")

kegg_host_fa_results <- perform_edger_and_plot(host_kegg_norm_imp, selected_matrix = "faeces", save = F,
                                               abundance_column = "intensity", save_name = "kegg_host")

table_p_host_fa_edger_kegg <- rbind(kegg_host_fa_results$`diet2-diet1`$table_sig, kegg_host_fa_results$`diet3-diet1`$table_sig, 
                                    kegg_host_fa_results$`diet4-diet1`$table_sig, kegg_host_fa_results$`diet3-diet2`$table_sig, 
                                    kegg_host_fa_results$`diet4-diet2`$table_sig, kegg_host_fa_results$`diet4-diet3`$table_sig)
write_tsv(table_p_host_fa_edger_kegg, "tables/93_table_p_host_fa_edger_kegg.txt")


## Multiomics
###

# load image
load("temp/73_multiomics_diablo_ileum.RData")
load("temp/74_multiomics_diablo_faeces.RData")
source(here("70_multiomics_functions.R"))

# associations for models in ileal digesta

pdf(NULL)
matrix_il_combined_1_2 <- circosPlot(diablo_il_combined_1_2, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

table_diablo_il_1_2 <- matrix_il_combined_1_2 %>%
  as_tibble(rownames = "feature")

write_tsv(table_diablo_il_1_2, "tables/93_table_diablo_il_1_2.txt")

pdf(NULL)
matrix_il_combined_1_3 <- circosPlot(diablo_il_combined_1_3, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

table_diablo_il_1_3 <- matrix_il_combined_1_3 %>%
  as_tibble(rownames = "feature")

write_tsv(table_diablo_il_1_3, "tables/93_table_diablo_il_1_3.txt")

pdf(NULL)
matrix_il_combined_1_4 <- circosPlot(diablo_il_combined_1_4, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

table_diablo_il_1_4 <- matrix_il_combined_1_4 %>%
  as_tibble(rownames = "feature")

write_tsv(table_diablo_il_1_4, "tables/93_table_diablo_il_1_4.txt")

# associations for models in feces

pdf(NULL)
matrix_fa_combined_1_2 <- circosPlot(diablo_fa_combined_1_2, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

table_diablo_fa_1_2 <- matrix_fa_combined_1_2 %>%
  as_tibble(rownames = "feature")

write_tsv(table_diablo_fa_1_2, "tables/93_table_diablo_fa_1_2.txt")

pdf(NULL)
matrix_fa_combined_1_3 <- circosPlot(diablo_fa_combined_1_3, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

table_diablo_fa_1_3 <- matrix_fa_combined_1_3 %>%
  as_tibble(rownames = "feature")

write_tsv(table_diablo_fa_1_3, "tables/93_table_diablo_fa_1_3.txt")

pdf(NULL)
matrix_fa_combined_1_4 <- circosPlot(diablo_fa_combined_1_4, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

table_diablo_fa_1_4 <- matrix_fa_combined_1_4 %>%
  as_tibble(rownames = "feature")

write_tsv(table_diablo_fa_1_4, "tables/93_table_diablo_fa_1_4.txt")