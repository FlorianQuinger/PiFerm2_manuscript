library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

#saveswitch
#save=TRUE
save =FALSE

# load meta

meta <- readRDS("clean/meta1.RDS") %>% filter(sampleid != "103")

# load proteins file with intensity

proteins_norm_imp <- readRDS("clean/4_proteins_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
proteins_norm_imp_rel_abd <- readRDS("clean/4_proteins_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

host_kegg_norm_imp <- readRDS("clean/4_host_kegg_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
host_kegg_norm_imp_rel_abd <- readRDS("clean/4_host_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

host_go_norm_imp <- readRDS("clean/4_host_go_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
host_go_norm_imp_rel_abd <- readRDS("clean/4_host_go_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

kegg_norm_imp <- readRDS("clean/4_kegg_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
kegg_norm_imp_rel_abd <- readRDS("clean/4_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

go_norm_imp <- readRDS("clean/4_go_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
go_norm_imp_rel_abd <- readRDS("clean/4_go_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

cog_norm_imp <- readRDS("clean/4_cog_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
cog_norm_imp_rel_abd <- readRDS("clean/4_cog_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

cazy_norm_imp <- readRDS("clean/4_cazy_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
cazy_norm_imp_rel_abd <- readRDS("clean/4_cazy_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

proteins_func_tax_combined <- readRDS("clean/4_function_taxa_combined_long_norm_imp_intensity_filtered.RDS") %>%
  inner_join(meta, by = "sampleid") %>% filter(sampleid != "103")

# load function files to merge
functions <- readRDS("clean/4_functions_long_raw.RDS") %>%
  mutate(cog_cat_name = str_c(cog_category, cog_name, sep = " - "))

host_functions <- readRDS("clean/4_host_functions_long_raw_intensity.RDS") %>%
  dplyr::select(-sampleid, -intensity) %>%
  distinct()

feed_functions <- readRDS("clean/4_feed_functions_long_raw_intensity.RDS") %>%
  dplyr::select(-sampleid, -intensity) %>%
  distinct()

# taxonomy file to assign proteins

taxonomy <- readRDS("clean/4_taxa_long_raw_intensity.RDS") %>%
  dplyr::select(-sampleid, -intensity) %>%
  distinct()

# split datasets

proteins_host <- proteins_norm_imp %>%
  filter(origin == "pig")
proteins_micro <- proteins_norm_imp %>%
  filter(origin == "micro")
proteins_feed <- proteins_norm_imp %>%
  filter(origin %in% c("soy", "pea", "rape", "barley", "wheat"))
proteins_pea <- proteins_norm_imp %>%
  filter(origin == "pea")

# try barplots

proteins_host_rel_abd <- proteins_norm_imp_rel_abd %>%
  filter(origin == "pig") %>%
  recalculate_rel_abd()
proteins_micro_rel_abd <- proteins_norm_imp_rel_abd %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()
proteins_feed_rel_abd <- proteins_norm_imp_rel_abd %>%
  filter(origin %in% c("soy", "pea", "rape", "barley", "wheat")) %>%
  recalculate_rel_abd()
proteins_pea_rel_abd <- proteins_norm_imp_rel_abd %>%
  filter(origin %in% c("pea")) %>%
  recalculate_rel_abd()

# for host proteins

taxa_barplot_from_ktable(dplyr::select(proteins_host_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ileal digesta", title = "Host proteins ileum",
                         save_name = "46_barplot_host_proteins_ileum", save = save)

taxa_barplot_from_ktable(dplyr::select(proteins_host_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "faeces", title = "Host proteins faeces",
                         save_name = "46_barplot_host_proteins_faeces", save = save, threshold = 1)

# for microbial proteins

taxa_barplot_from_ktable(dplyr::select(proteins_micro_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ileal digesta", title = "Microbial proteins ileum",
                         save_name = "46_barplot_micro_proteins_ileum", save = save, threshold = 0.4)

taxa_barplot_from_ktable(dplyr::select(proteins_micro_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "faeces", title = "Microbial proteins faeces",
                         save_name = "46_barplot_micro_proteins_faeces", save = save, threshold = 0.3)

# for feed proteins

taxa_barplot_from_ktable(dplyr::select(proteins_feed_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ileal digesta", title = "Feed proteins ileum",
                         save_name = "46_barplot_feed_proteins_ileum", save = save, threshold = 2)

taxa_barplot_from_ktable(dplyr::select(proteins_feed_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "faeces", title = "Feed proteins faeces",
                         save_name = "46_barplot_feed_proteins_faeces", save = save, threshold = 2)

# for pea proteins

taxa_barplot_from_ktable(dplyr::select(proteins_pea_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "ileal digesta", title = "Pea proteins ileum",
                         save_name = "46_barplot_pea_proteins_ileum", save = save, threshold = 2)

taxa_barplot_from_ktable(dplyr::select(proteins_pea_rel_abd, name = proteinid, sampleid, rel_abd),meta = meta,
                         selected_rank = "low", selected_matrix = "faeces", title = "Pea proteins faeces",
                         save_name = "46_barplot_pea_proteins_faeces", save = save, threshold = 2)

# perform analysis


proteins_micro_il_results <- perform_edger_and_plot(proteins_micro, selected_matrix = "ileum", save = T,
                                              abundance_column = "intensity", save_name = "proteingroups_micro")

proteins_micro_fa_results <- perform_edger_and_plot(proteins_micro, selected_matrix = "faeces", save = T,
                                              abundance_column = "intensity", save_name = "proteingroups_micro")

proteins_host_il_results <- perform_edger_and_plot(proteins_host, selected_matrix = "ileum", save = T,
                                              abundance_column = "intensity", save_name = "proteingroups_host")

proteins_host_fa_results <- perform_edger_and_plot(proteins_host, selected_matrix = "faeces", save = T,
                                              abundance_column = "intensity", save_name = "proteingroups_host")

proteins_feed_il_results <- perform_edger_and_plot(proteins_feed, selected_matrix = "ileum", save = T,
                                                   abundance_column = "intensity", save_name = "proteingroups_feed")

proteins_feed_fa_results <- perform_edger_and_plot(proteins_feed, selected_matrix = "faeces", save = T,
                                                   abundance_column = "intensity", save_name = "proteingroups_feed")

proteins_pea_il_results <- perform_edger_and_plot(proteins_pea, selected_matrix = "ileum", save = T,
                                                   abundance_column = "intensity", save_name = "proteingroups_feed")

proteins_pea_fa_results <- perform_edger_and_plot(proteins_pea, selected_matrix = "faeces", save = T,
                                                   abundance_column = "intensity", save_name = "proteingroups_feed")

kegg_il_results <- perform_edger_and_plot(kegg_norm_imp, selected_matrix = "ileum", save = T,
                                                   abundance_column = "intensity", save_name = "kegg_micro")

kegg_fa_results <- perform_edger_and_plot(kegg_norm_imp, selected_matrix = "faeces", save = T,
                                                   abundance_column = "intensity", save_name = "kegg_micro")

go_il_results <- perform_edger_and_plot(go_norm_imp, selected_matrix = "ileum", save = T,
                                          abundance_column = "intensity", save_name = "go_micro")

go_fa_results <- perform_edger_and_plot(go_norm_imp, selected_matrix = "faeces", save = T,
                                          abundance_column = "intensity", save_name = "go_micro")

cog_il_results <- perform_edger_and_plot(cog_norm_imp, selected_matrix = "ileum", save = T,
                                        abundance_column = "intensity", save_name = "cog_micro")

cog_fa_results <- perform_edger_and_plot(cog_norm_imp, selected_matrix = "faeces", save = T,
                                        abundance_column = "intensity", save_name = "cog_micro")

cazy_il_results <- perform_edger_and_plot(cazy_norm_imp, selected_matrix = "ileum", save = T,
                                         abundance_column = "intensity", save_name = "cazy_micro") # no sig at all

cazy_fa_results <- perform_edger_and_plot(cazy_norm_imp, selected_matrix = "faeces", save = T,
                                         abundance_column = "intensity", save_name = "cazy_micro")

kegg_host_il_results <- perform_edger_and_plot(host_kegg_norm_imp, selected_matrix = "ileum", save = T,
                                          abundance_column = "intensity", save_name = "kegg_host")

kegg_host_fa_results <- perform_edger_and_plot(host_kegg_norm_imp, selected_matrix = "faeces", save = T,
                                          abundance_column = "intensity", save_name = "kegg_host")

go_host_il_results <- perform_edger_and_plot(host_go_norm_imp, selected_matrix = "ileum", save = T,
                                        abundance_column = "intensity", save_name = "go_host")

go_host_fa_results <- perform_edger_and_plot(host_go_norm_imp, selected_matrix = "faeces", save = T,
                                        abundance_column = "intensity", save_name = "go_host")

#add taxonomy to micro proteins
proteins_micro_il_results <- add_taxonomy(proteins_micro_il_results, taxonomy, join_by = "bin")
proteins_micro_fa_results <- add_taxonomy(proteins_micro_fa_results, taxonomy, join_by = "bin")

# add functional annotation to result lists

# add to host proteins
proteins_host_il_results <- add_functions(proteins_host_il_results, host_functions, join_by = "proteinid")
proteins_host_fa_results <- add_functions(proteins_host_fa_results, host_functions, join_by = "proteinid")

#add to micro proteins
proteins_micro_il_results <- add_functions(proteins_micro_il_results, functions, join_by = "proteinid")
proteins_micro_fa_results <- add_functions(proteins_micro_fa_results, functions, join_by = "proteinid")

# add to feed proteins 

proteins_feed_il_results <- add_functions(proteins_feed_il_results, feed_functions, join_by = "proteinid")
proteins_feed_fa_results <- add_functions(proteins_feed_fa_results, feed_functions, join_by = "proteinid")

# add to pea proteins

proteins_pea_il_results <- add_functions(proteins_pea_il_results, feed_functions, join_by = "proteinid")
proteins_pea_fa_results <- add_functions(proteins_pea_fa_results, feed_functions, join_by = "proteinid")

# KEGG

enrich_kegg(kegg_il_results, save_name = "46_enrich_kegg_ileum", save = save)
enrich_kegg(kegg_fa_results, save_name = "46_enrich_kegg_faeces", save = save)

enrich_kegg(proteins_micro_il_results_kegg, save_name = "46_enrich_kegg_proteingroups_ileum", save = save)
enrich_kegg(proteins_micro_fa_results_kegg, save_name = "46_enrich_kegg_proteingroups_faeces", save = save)

# kegg module enrichment

enrich_mkegg(proteins_micro_il_results_kegg, save_name = "46_enrich_mkegg_proteingroups_ileum", save = save)
enrich_mkegg(proteins_micro_fa_results_kegg, save_name = "46_enrich_mkegg_proteingroups_faeces", save = save)

enrich_kegg(kegg_il_results, save_name = "46_enrich_mkegg_ileum", save = save)
enrich_kegg(kegg_fa_results, save_name = "46_enrich_mkegg_faeces", save = save)

#COG

enrich_cog(cog_il_results, save_name = "46_enrichment_cog_ileum") 
enrich_cog(cog_fa_results, save_name = "46_enrichment_cog_faeces")

enrich_cog(proteins_micro_il_results_cog, save_name = "46_enrichment_cog_proteingroups_ileum", save = save)
enrich_cog(proteins_micro_fa_results_cog, save_name = "46_enrichment_cog_proteingroups_faeces", save = save)