library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

#save switch
#save = TRUE
save = FALSE

# load meta

meta <- readRDS("clean/meta1.RDS") %>% filter(sampleid != "103")

# load files

k2report_rel_abd <- readRDS("clean/3_k2report_reads_long_rel_abd.RDS")%>% filter(sampleid != "103")
k2report_reads <- readRDS("clean/3_k2report_reads_long_reads.RDS")%>% filter(sampleid != "103")
breport_rel_abd <- readRDS("clean/3_breport_reads_long_rel_abd.RDS")%>% filter(sampleid != "103")
breport_reads <- readRDS("clean/3_breport_reads_long_reads.RDS")%>% filter(sampleid != "103")
bracken_rel_abd <- readRDS("clean/3_bracken_reads_long_rel_abd.RDS")%>% filter(sampleid != "103")
bracken_reads <- readRDS("clean/3_bracken_reads_long_reads.RDS")%>% filter(sampleid != "103")

# load filtered files

breport_rel_abd_filtered <- readRDS("clean/3_breport_reads_long_rel_abd_filtered.RDS")%>% filter(sampleid != "103")
breport_reads_filtered <- readRDS("clean/3_breport_reads_long_reads_filtered.RDS")%>% filter(sampleid != "103")
bracken_rel_abd_filtered <- readRDS("clean/3_bracken_reads_long_rel_abd_filtered.RDS")%>% filter(sampleid != "103")
bracken_reads_filtered <- readRDS("clean/3_bracken_reads_long_reads_filtered.RDS")%>% filter(sampleid != "103")

#########

breport_reads %>% filter(sampleid == "101", rank == "R1") %>% head() # bracken does not assign all reads, but does not show the number of unclassified reads
k2report_reads %>% filter(sampleid == "101", rank == "R1") %>% head() 

# Unclassified reads by Kraken

k2report_rel_abd_unclassified <- k2report_reads %>%
  dplyr::select(-rank) %>%
  filter(name %in% c("root", "unclassified")) %>%
  pivot_wider(names_from = name, values_from = reads) %>%
  mutate(rel_abd = unclassified / (unclassified + root)) %>% # total is unclassified+root
  inner_join(dplyr::select(meta, sampleid, matrix), by = "sampleid")

k2report_rel_abd_unclassified %>%
  group_by(matrix) %>%
  summarise(mean(rel_abd))

ggplot(k2report_rel_abd_unclassified, aes(x = matrix, y = rel_abd)) +
  geom_violin(fill = "lightblue") +
  geom_quasirandom(size = 3, width = 0.3) +
  labs(x = "", y = "reads", title = "unclassified reads by Kraken2")

save_big("32_kraken2_unclassified")

# Alpha diversity

## observed species

observed_species <- bracken_reads %>%
  mutate(reads = ifelse(reads > 0, 1, 0)) %>% # set to absence presence
  group_by(sampleid) %>%
  summarise(observed = sum(reads))

## shannon

bracken_matrix <- bracken_rel_abd %>%
  pivot_wider(names_from = "name", values_from = "rel_abd") %>%
  column_to_rownames("sampleid") 

shannon <- tibble(sampleid = rownames(bracken_matrix),
                  shannon = diversity(bracken_matrix, index = "shannon"),
                  simpson = diversity(bracken_matrix, index = "invsimpson"))

# plot

alpha <- observed_species %>%
  inner_join(shannon, by = "sampleid") %>%
  inner_join(meta, by = "sampleid")

ggplot(pivot_longer(alpha, c(shannon, observed, simpson), names_to = "index", values_to = "value"), aes(x = diet, y = value)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom() +
  facet_wrap(vars(matrix, index), scales = "free_y")

# comparison 

alpha_il <- select_matrix(alpha, "ileal digesta")
alpha_fa <- select_matrix(alpha, "faeces")

comp_obs_il <- combined_comparison(select_response(alpha_il, "observed"), transformation = "test")
comp_shannon_il <- combined_comparison(select_response(alpha_il, "shannon"), transformation = "test")
comp_simpson_il <- combined_comparison(select_response(alpha_il, "simpson"), transformation = "test")
comp_obs_fa <- combined_comparison(select_response(alpha_fa, "observed"), transformation = "test")
comp_shannon_fa <- combined_comparison(select_response(alpha_fa, "shannon"), transformation = "test")
comp_simpson_fa <- combined_comparison(select_response(alpha_fa, "simpson"), transformation = "test")

table_obs_il <- create_results_table(select_response(alpha_il, "observed"), comp_obs_il, response = "observed_il")
table_shannon_il <- create_results_table(select_response(alpha_il, "shannon"), comp_shannon_il, response = "shannon_il")
table_simpson_il <- create_results_table(select_response(alpha_il, "simpson"), comp_simpson_il, response = "simpson_il")
table_obs_fa <- create_results_table(select_response(alpha_fa, "observed"), comp_obs_fa, response = "observed_fa")
table_shannon_fa <- create_results_table(select_response(alpha_fa, "shannon"), comp_shannon_fa, response = "shannon_fa")
table_simpson_fa <- create_results_table(select_response(alpha_fa, "simpson"), comp_simpson_fa, response = "simpson_fa")

table_alpha <- table_obs_il %>%
  inner_join(table_shannon_il, by = "diet") %>%
  inner_join(table_simpson_il, by = "diet") %>%
  inner_join(table_obs_fa, by = "diet") %>%
  inner_join(table_shannon_fa, by = "diet") %>%
  inner_join(table_simpson_fa, by = "diet")

write_tsv(table_alpha, "tables/32_reads_alpha_div.txt")


# Beta diversity with unfiltered frame

do_ordination(bracken_rel_abd)

do_ordination(bracken_rel_abd, region = "ileum", title = "Taxonomy reads ileum", save_name = "ileum_reads", save = save, second_indicator = "animal")
#do_ordination(bracken_rel_abd, region = "ileum") # subtile changes by filtering

do_ordination(bracken_rel_abd, region = "faeces", title = "Taxonomy reads faeces", save_name = "faeces_reads", save = save, second_indicator = "animal")
#do_ordination(bracken_rel_abd, region = "faeces") # subtile changes by filtering

# Taxa barplots ileum

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "P", selected_matrix = "ileal digesta", 
                         title = "Phyla ileal digesta", save_name = "32_taxa_barplot_ileum_reads_phylum", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "C", selected_matrix = "ileal digesta", 
                         title = "Class ileal digesta", save_name = "32_taxa_barplot_ileum_reads_class", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "O", selected_matrix = "ileal digesta", 
                         title = "Orders ileal digesta", save_name = "32_taxa_barplot_ileum_reads_order", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "F", selected_matrix = "ileal digesta", 
                         title = "Families ileal digesta", save_name = "32_taxa_barplot_ileum_reads_family", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "G", selected_matrix = "ileal digesta", 
                         title = "Genera ileal digesta", save_name = "32_taxa_barplot_ileum_reads_genus", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "ileal digesta", 
                         title = "Species ileal digesta", save_name = "32_taxa_barplot_ileum_reads_species", save = save)

breport_rel_abd_filtered %>%
  filter_ileum() %>%
  filter(rank == "G") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  arrange(desc(rel_abd))

breport_rel_abd_filtered %>%
  filter_faeces() %>%
  filter(rank == "G") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  arrange(desc(rel_abd))

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "ileal digesta", 
                         title = "Species ileal digesta", save_name = "32_taxa_barplot_ileum_reads_species_animal_effect", save = save,
                         grouping_factors = c("animal", "diet"))

#taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "ileal digesta") # slight changes for filtering

# Taxa barplots faeces

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "P", selected_matrix = "faeces", 
                         title = "Phyla faeces", save_name = "32_taxa_barplot_faeces_reads_phylum", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "C", selected_matrix = "faeces", 
                         title = "Class faeces", save_name = "32_taxa_barplot_faeces_reads_class", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "O", selected_matrix = "faeces", 
                         title = "Orders faeces", save_name = "32_taxa_barplot_faeces_reads_order", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "F", selected_matrix = "faeces", 
                         title = "Families faeces", save_name = "32_taxa_barplot_faeces_reads_family", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "G", selected_matrix = "faeces", 
                         title = "Genera faeces", save_name = "32_taxa_barplot_faeces_reads_genus", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "faeces", threshold = 0.5, 
                         title = "Species faeces", save_name = "32_taxa_barplot_faeces_reads_species", save = save)

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "faeces", threshold = .75, 
                         title = "Species faeces", save_name = "32_taxa_barplot_faeces_reads_species_animal_effect", save = save,
                         grouping_factors = c("animal", "diet"))

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "faeces", threshold = .75, 
                         title = "Species faeces", save = save,
                         grouping_factors = c("diet", "animal"))

#taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "S", selected_matrix = "faeces", threshold = 0.5)

# differential abundance

# ileal digesta

il_ancom_P <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "P", 
                                       selected_matrix = "ileal digesta")
il_ancom_C <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "C", 
                                       selected_matrix = "ileal digesta")
il_ancom_O <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "O", 
                                       selected_matrix = "ileal digesta")
il_ancom_F <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "F", 
                                       selected_matrix = "ileal digesta")
il_ancom_G <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "G", 
                                       selected_matrix = "ileal digesta")
il_ancom_S <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "S", 
                                       selected_matrix = "ileal digesta")

il_ancom_low <- perform_ancombc_and_plot(ktable = bracken_reads_filtered, meta = meta, selected_rank = "low", 
                                       selected_matrix = "ileal digesta")

# faeces

fa_ancom_P <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "P", 
                                       selected_matrix = "faeces")
fa_ancom_C <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "C", 
                                       selected_matrix = "faeces")
fa_ancom_O <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "O", 
                                       selected_matrix = "faeces")
fa_ancom_F <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "F", 
                                       selected_matrix = "faeces")
fa_ancom_G <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "G", 
                                       selected_matrix = "faeces")
fa_ancom_S <- perform_ancombc_and_plot(ktable = breport_reads_filtered, meta = meta, selected_rank = "S", 
                                       selected_matrix = "faeces")

fa_ancom_low <- perform_ancombc_and_plot(ktable = bracken_reads_filtered, meta = meta, selected_rank = "low", 
                                         selected_matrix = "faeces")

# save significant taxa for correlation

il_p_sig <- create_venn_for_sig(il_ancom_P, breport_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "P", save_name = "tax", save = save)
il_c_sig <- create_venn_for_sig(il_ancom_C, breport_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "C", save_name = "tax", save = save)
il_o_sig <- create_venn_for_sig(il_ancom_O, breport_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "O", save_name = "tax", save = save)
il_f_sig <- create_venn_for_sig(il_ancom_F, breport_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "F", save_name = "tax", save = save)
il_g_sig <- create_venn_for_sig(il_ancom_G, breport_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "G", save_name = "tax", save = save)
il_s_sig <- create_venn_for_sig(il_ancom_S, breport_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "S", save_name = "tax", save = save)
il_low_sig <- create_venn_for_sig(il_ancom_low, bracken_rel_abd_filtered, selected_matrix = "ileal digesta",
                                selected_rank = "low", save_name = "tax", save = save)

fa_p_sig <- create_venn_for_sig(fa_ancom_P, breport_rel_abd_filtered, selected_matrix = "faeces",
                                selected_rank = "P", save_name = "tax", save = save)
fa_c_sig <- create_venn_for_sig(fa_ancom_C, breport_rel_abd_filtered, selected_matrix = "faeces",
                                selected_rank = "C", save_name = "tax", save = save)
fa_o_sig <- create_venn_for_sig(fa_ancom_O, breport_rel_abd_filtered, selected_matrix = "faeces",
                                selected_rank = "O", save_name = "tax", save = save)
fa_f_sig <- create_venn_for_sig(fa_ancom_F, breport_rel_abd_filtered, selected_matrix = "faeces",
                                selected_rank = "F", save_name = "tax", save = save)
fa_g_sig <- create_venn_for_sig(fa_ancom_G, breport_rel_abd_filtered, selected_matrix = "faeces",
                                selected_rank = "G", save_name = "tax", save = save)
fa_s_sig <- create_venn_for_sig(fa_ancom_S, breport_rel_abd_filtered, selected_matrix = "faeces",
                                selected_rank = "S", save_name = "tax", save = save)
fa_low_sig <- create_venn_for_sig(fa_ancom_low, bracken_rel_abd_filtered, selected_matrix = "faeces",
                                  selected_rank = "low", save_name = "tax", save = save)

# create plots for sig taxa

il_s_sig_median <- il_s_sig %>%
  group_by(name) %>%
  summarise(rel_abd = median(rel_abd))


fa_s_sig_median <- fa_s_sig %>%
  group_by(name) %>%
  summarise(rel_abd = median(rel_abd))

#ileal digesta

create_plot_for_taxon(output_object = il_ancom_F, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "F", 
                      selected_matrix = "ileal digesta", selected_taxon = "Campylobacteraceae")

create_plot_for_taxon(output_object = il_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "ileal digesta", selected_taxon = "Prevotella pectinovora")

create_plot_for_taxon(output_object = il_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "ileal digesta", selected_taxon = "Lachnospira hominis")

create_plot_for_taxon(output_object = il_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "ileal digesta", selected_taxon = "Selenomonas_C bovis")


#faeces

create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "faeces", selected_taxon = "Bifidobacterium boum")

create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "faeces", selected_taxon = "Cryptobacteroides sp900546925")

create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "faeces", selected_taxon = "Fusicatenibacter saccharivorans")

create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "faeces", selected_taxon = "Methanobrevibacter_A sp022775905")

create_plot_for_taxon(output_object = fa_ancom_G, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "G", 
                      selected_matrix = "faeces", selected_taxon = "UBA2810")

create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "faeces", selected_taxon = "UBA2810 sp002351705")

create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered, selected_rank = "S", 
                      selected_matrix = "faeces", selected_taxon = "UBA2810 sp900317945")

# loop over significant genera in ileum and faeces and plot them

for (i in 1:length(unique(il_s_sig$name))) {
  taxa <- unique(il_s_sig$name)[i]
  create_plot_for_taxon(output_object = il_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered,
                        selected_rank = "S", selected_matrix = "ileal digesta", selected_taxon = taxa)
}


for (i in 1:length(unique(fa_f_sig$name))) {
  taxa <- unique(fa_f_sig$name)[i]
  create_plot_for_taxon(output_object = fa_ancom_F, meta = meta, df_rel_abd = breport_rel_abd_filtered,
                        selected_rank = "F", selected_matrix = "faeces", selected_taxon = taxa)
}

for (i in 1:length(unique(fa_g_sig$name))) {
  taxa <- unique(fa_g_sig$name)[i]
  create_plot_for_taxon(output_object = fa_ancom_G, meta = meta, df_rel_abd = breport_rel_abd_filtered,
                        selected_rank = "G", selected_matrix = "faeces", selected_taxon = taxa)
}

for (i in 1:length(unique(fa_s_sig$name))) {
  taxa <- unique(fa_s_sig$name)[i]
  create_plot_for_taxon(output_object = fa_ancom_S, meta = meta, df_rel_abd = breport_rel_abd_filtered,
                        selected_rank = "S", selected_matrix = "faeces", selected_taxon = taxa)
}


