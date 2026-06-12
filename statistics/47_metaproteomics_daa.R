library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta

meta <- readRDS("clean/meta1.RDS") %>% filter(sampleid != "103")

#saveswitch
#save=TRUE
save =FALSE

# load taxonomy files

taxonomy_norm_imp_intensity <- readRDS("clean/4_taxonomy_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
taxa_norm_imp_intensity <- readRDS("clean/4_taxa_long_norm_imp_intensity_filtered.RDS") %>%
  dplyr::select(bin, sampleid, intensity)  %>% filter(sampleid != "103")
taxa_raw_rel_abd <- readRDS("clean/4_taxa_long_raw_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

taxonomy_norm_imp_rel_abd <- readRDS("clean/4_taxonomy_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
taxa_norm_imp_rel_abd <- readRDS("clean/4_taxa_long_norm_imp_rel_abd_filtered.RDS") %>%
  dplyr::select(bin, sampleid, rel_abd)  %>% filter(sampleid != "103")

# Alpha diversity

## observed species

observed_species <- taxa_raw_rel_abd %>%
  mutate(rel_abd = ifelse(rel_abd > 0, 1, 0)) %>% # set to absence presence
  group_by(sampleid) %>%
  summarise(observed = sum(rel_abd))

## shannon

matrix <- taxa_raw_rel_abd %>%
  dplyr::select(bin, sampleid, rel_abd) %>%
  pivot_wider(names_from = "bin", values_from = "rel_abd") %>%
  column_to_rownames("sampleid") %>%
  mutate_all(~ifelse(is.na(.), 0, .))

shannon <- tibble(sampleid = rownames(matrix),
                  shannon = diversity(matrix, index = "shannon"),
                  simpson = diversity(matrix, index = "invsimpson"))

# plot

alpha <- observed_species %>%
  inner_join(shannon, by = "sampleid") %>%
  inner_join(meta, by = "sampleid")

ggplot(pivot_longer(alpha, c(shannon, observed, simpson), names_to = "index", values_to = "value"), aes(x = diet, y = value)) +
  geom_boxplot(outliers = F) +
  geom_quasirandom() +
  facet_wrap(vars(matrix, index), scales = "free_y")

#comparison

alpha_il <- select_matrix(alpha, "ileal digesta")
alpha_fa <- select_matrix(alpha, "faeces")

comp_shannon_il <- combined_comparison(select_response(alpha_il, "shannon"), transformation = "test") # not sig
comp_obs_il <- combined_comparison(select_response(alpha_il, "observed"), transformation = "test")
comp_simpson_il <- combined_comparison(select_response(alpha_il, "simpson"), transformation = "test")
comp_shannon_fa <- combined_comparison(select_response(alpha_fa, "shannon"), transformation = "test") #not sig
comp_simpson_fa <- combined_comparison(select_response(alpha_fa, "simpson"), transformation = "test")

table_obs_il <- create_results_table(select_response(alpha_il, "observed"), comp_obs_il, response = "observed_il")
table_shannon_il <- create_results_table(select_response(alpha_il, "shannon"), comp_shannon_il, response = "shannon_il")
table_simpson_il <- create_results_table(select_response(alpha_il, "simpson"), comp_simpson_il, response = "simpson_il")
table_shannon_fa <- create_results_table(select_response(alpha_fa, "shannon"), comp_shannon_fa, response = "shannon_fa")
table_simpson_fa <- create_results_table(select_response(alpha_fa, "simpson"), comp_simpson_fa, response = "simpson_fa")

table_alpha <- table_obs_il %>%
  inner_join(table_shannon_il, by = "diet") %>%
  inner_join(table_simpson_il, by = "diet") %>%
  #inner_join(table_obs_fa, by = "diet") %>%
  inner_join(table_shannon_fa, by = "diet") %>%
  inner_join(table_simpson_fa, by = "diet")

write_tsv(table_alpha, "tables/47_alpha_diversity.txt")

# Taxa barplots

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "P", selected_matrix = "ileal digesta", 
                         title = "Phyla ileal digesta", save_name = "47_taxa_barplot_ileum_phylum", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "C", selected_matrix = "ileal digesta", 
                         title = "Class ileal digesta", save_name = "47_taxa_barplot_ileum_class", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "O", selected_matrix = "ileal digesta", 
                         title = "Orders ileal digesta", save_name = "47_taxa_barplot_ileum_order", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "F", selected_matrix = "ileal digesta", 
                         title = "Family ileal digesta", save_name = "47_taxa_barplot_ileum_family", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", selected_matrix = "ileal digesta", 
                         title = "Genus ileal digesta", save_name = "47_taxa_barplot_ileum_genus", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", selected_matrix = "ileal digesta", grouping_factors = c("animal", "diet"),
                         title = "Genus ileal digesta", save_name = "47_taxa_barplot_ileum_genus_animal_effect", save = save)

taxonomy_norm_imp_rel_abd %>%
  filter_ileum() %>%
  filter(rank == "G") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  arrange(desc(rel_abd))

taxonomy_norm_imp_rel_abd %>%
  filter_faeces() %>%
  filter(rank == "G") %>%
  group_by(name) %>%
  summarise(rel_abd = mean(rel_abd)) %>%
  arrange(desc(rel_abd))
  
taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "S", selected_matrix = "ileal digesta", 
                         title = "Species ileal digesta", save_name = "47_taxa_barplot_ileum_species", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "S", selected_matrix = "ileal digesta", grouping_factors = c("animal", "diet"), threshold = 1.5,
                         title = "Species ileal digesta", save_name = "47_taxa_barplot_ileum_species_animal_effect", save = save)

# Faeces

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "P", selected_matrix = "faeces", 
                         title = "Phyla faeces", save_name = "47_taxa_barplot_faeces_phylum", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "C", selected_matrix = "faeces", 
                         title = "Class faeces", save_name = "47_taxa_barplot_faeces_class", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "O", selected_matrix = "faeces", 
                         title = "Orders faeces", save_name = "47_taxa_barplot_faeces_order", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "F", selected_matrix = "faeces", 
                         title = "Family faeces", save_name = "47_taxa_barplot_faeces_family", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", selected_matrix = "faeces", 
                         title = "Genus faeces", save_name = "47_taxa_barplot_faeces_genus", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", selected_matrix = "faeces", 
                         grouping_factors = c("animal", "diet"),
                         title = "Genus faeces", save_name = "47_taxa_barplot_faeces_genus_animal_effect", save = save)

taxa_barplot_from_ktable(taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "S", selected_matrix = "faeces", 
                         title = "Species faeces", save_name = "47_taxa_barplot_faeces_species", save = save)

# differential abundance analysis

# ileum

il_ancom_P <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "P", 
                                       selected_matrix = "ileal digesta")

il_ancom_C <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "C", 
                                       selected_matrix = "ileal digesta")

il_ancom_O <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "O", 
                                       selected_matrix = "ileal digesta")

il_ancom_F <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "F", 
                                       selected_matrix = "ileal digesta")

il_ancom_G <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "G", 
                                       selected_matrix = "ileal digesta")

il_ancom_S <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "S", 
                                       selected_matrix = "ileal digesta")

il_ancom_low <- perform_ancombc_and_plot(ktable = taxa_norm_imp_intensity, meta = meta, selected_rank = "low", 
                                       selected_matrix = "ileal digesta")

# plot

create_plot_for_taxon(il_ancom_S, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "S",
                      selected_matrix = "ileal digesta", selected_taxon = "UBA4248 sp945875805")


create_plot_for_taxon(il_ancom_S, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "S",
                      selected_matrix = "ileal digesta", selected_taxon = "Prevotella sp002299635")

# faeces

fa_ancom_P <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "P", 
                                       selected_matrix = "faeces")

fa_ancom_C <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "C", 
                                       selected_matrix = "faeces")

fa_ancom_O <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "O", 
                                       selected_matrix = "faeces")

fa_ancom_F <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "F", 
                                       selected_matrix = "faeces")

fa_ancom_G <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "G", 
                                       selected_matrix = "faeces")

fa_ancom_S <- perform_ancombc_and_plot(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "S", 
                                       selected_matrix = "faeces")

fa_ancom_low <- perform_ancombc_and_plot(ktable = taxa_norm_imp_intensity, meta = meta, selected_rank = "low", 
                                         selected_matrix = "faeces")

# save significant 

il_p_sig <- create_venn_for_sig(il_ancom_P, taxonomy_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                selected_rank = "P", save_name = "tax", save = save)
il_c_sig <- create_venn_for_sig(il_ancom_C, taxonomy_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                selected_rank = "C", save_name = "tax", save = save)
il_o_sig <- create_venn_for_sig(il_ancom_O, taxonomy_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                selected_rank = "O", save_name = "tax", save = save)
il_f_sig <- create_venn_for_sig(il_ancom_F, taxonomy_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                selected_rank = "F", save_name = "tax", save = save)
il_g_sig <- create_venn_for_sig(il_ancom_G, taxonomy_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                selected_rank = "G", save_name = "tax", save = save)
il_s_sig <- create_venn_for_sig(il_ancom_S, taxonomy_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                selected_rank = "S", save_name = "tax", save = save)
il_low_sig <- create_venn_for_sig(il_ancom_low, taxa_norm_imp_rel_abd, selected_matrix = "ileal digesta",
                                  selected_rank = "low", save_name = "tax", save = save)

fa_p_sig <- create_venn_for_sig(fa_ancom_P, taxonomy_norm_imp_rel_abd, selected_matrix = "faeces",
                                selected_rank = "P", save_name = "tax", save = save)
fa_c_sig <- create_venn_for_sig(fa_ancom_C, taxonomy_norm_imp_rel_abd, selected_matrix = "faeces",
                                selected_rank = "C", save_name = "tax", save = save)
fa_o_sig <- create_venn_for_sig(fa_ancom_O, taxonomy_norm_imp_rel_abd, selected_matrix = "faeces",
                                selected_rank = "O", save_name = "tax", save = save)
fa_f_sig <- create_venn_for_sig(fa_ancom_F, taxonomy_norm_imp_rel_abd, selected_matrix = "faeces",
                                selected_rank = "F", save_name = "tax", save = save)
fa_g_sig <- create_venn_for_sig(fa_ancom_G, taxonomy_norm_imp_rel_abd, selected_matrix = "faeces",
                                selected_rank = "G", save_name = "tax", save = save)
fa_s_sig <- create_venn_for_sig(fa_ancom_S, taxonomy_norm_imp_rel_abd, selected_matrix = "faeces",
                                selected_rank = "S", save_name = "tax", save = save)
fa_low_sig <- create_venn_for_sig(fa_ancom_low, taxa_norm_imp_rel_abd, selected_matrix = "faeces",
                                  selected_rank = "low", save_name = "tax", save = save)

# create individual plots for taxa

fa_g_sig_median <- fa_g_sig %>%
  group_by(name) %>%
  summarise(rel_abd = median(rel_abd))

create_plot_for_taxon(fa_ancom_G, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G",
                      selected_matrix = "faeces", selected_taxon = "Selenomonas_C", save = save)

create_plot_for_taxon(fa_ancom_G, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G",
                      selected_matrix = "faeces", selected_taxon = "Cryptobacteroides", save = save)

create_plot_for_taxon(fa_ancom_G, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G",
                      selected_matrix = "faeces", selected_taxon = "Ruminococcus", save = save)

create_plot_for_taxon(fa_ancom_G, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G",
                      selected_matrix = "faeces", selected_taxon = "Escherichia", save = save)

create_plot_for_taxon(fa_ancom_G, df_rel_abd = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G",
                      selected_matrix = "faeces", selected_taxon = "Bruticola", save = save)

for (i in 1:length(unique(fa_g_sig$name))) {
  taxa <- unique(fa_g_sig$name)[i]
  create_plot_for_taxon(output_object = fa_ancom_G, meta = meta, df_rel_abd = taxonomy_norm_imp_rel_abd,
                        selected_rank = "G", selected_matrix = "faeces", selected_taxon = taxa)
}

for (i in 1:length(unique(il_s_sig$name))) {
  taxa <- unique(il_s_sig$name)[i]
  create_plot_for_taxon(output_object = il_ancom_S, meta = meta, df_rel_abd = taxonomy_norm_imp_rel_abd,
                        selected_rank = "S", selected_matrix = "ileal digesta", selected_taxon = taxa)
}

