library(here)
library(venn)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta
dry_matter <- read_tsv("data/dry_matter.txt") %>%
  mutate(sampleid = as.character(sampleid))

meta <- readRDS("clean/meta1.RDS") %>%
  filter(sampleid != "103") %>%
  inner_join(dry_matter, by = "sampleid")

# load nmr files

nmr_ileum_dm <- readRDS("clean/5_nmr_ileum_wide_dm.RDS") %>%
  filter(sampleid != "103")
nmr_feces_dm <- readRDS("clean/5_nmr_feces_wide_dm.RDS")
nmr_ileum_fm <- readRDS("clean/5_nmr_ileum_wide_fm.RDS") %>%
  filter(sampleid != "103")
nmr_feces_fm <- readRDS("clean/5_nmr_feces_wide_fm.RDS")

# saveswitch

save = FALSE
# save = TRUE


# prepare dfs for analysis function

il_dm <- prepare_metabolomics(nmr_ileum_dm)
il_fm <- prepare_metabolomics(nmr_ileum_fm)
fe_dm <- prepare_metabolomics(nmr_feces_dm)
fe_fm <- prepare_metabolomics(nmr_feces_fm)


# ph
il_ph <- meta %>%
  filter(matrix == "ileal digesta") %>%
  dplyr::select(sampleid, response = ph) %>%
  prepare_metabolomics()

il_ph_out <- combined_comparison(il_ph, transformation = "test") # no effect on ph
plot_pairwise(filtered_df = il_ph, output = il_ph_out, selected_response = "ph")
save_big("52_ileum_ph")
il_ph_table <- create_results_table(il_ph, il_ph_out, response_name = "pH ileum")

fe_ph <- meta %>%
  filter(matrix == "faeces") %>%
  dplyr::select(sampleid, response = ph) %>%
  prepare_metabolomics()

fe_ph_out <- combined_comparison(fe_ph, transformation = "test") # effect of treatment on ph -> batch effects of etoh
create_results_plot(fe_ph, fe_ph_out, "pH feces", y_axis = "pH")
save_big("52_result_faeces_ph")
fe_ph_table <- create_results_table(fe_ph, fe_ph_out, response_name = "pH feces")

ph_table <- inner_join(il_ph_table, fe_ph_table, by = "diet")

# DM
il_dm <- meta %>%
  filter(matrix == "ileal digesta") %>%
  dplyr::select(sampleid, response = dry_matter) %>%
  prepare_metabolomics()

il_dm_out <- combined_comparison(il_dm, transformation = "test") # no effect
plot_pairwise(filtered_df = il_dm, output = il_dm_out, selected_response = "dry_matter")
save_big("52_ileum_dm")
il_dm_table <- create_results_table(il_dm, il_dm_out, response_name = "DM ileum")

fe_dm <- meta %>%
  filter(matrix == "faeces") %>%
  dplyr::select(sampleid, response = "dry_matter") %>%
  prepare_metabolomics()

fe_dm_out <- combined_comparison(fe_dm, transformation = "test") # no effect on DM
create_results_plot(fe_dm, fe_dm_out, "DM faeces", y_axis = "dry_matter")
save_big("52_result_faeces_dm")
fe_dm_table <- create_results_table(fe_dm, fe_dm_out, response_name = "DM feces")

dm_table <- inner_join(il_dm_table, fe_dm_table, by = "diet")


# dry matter

il_dm_table <- loop_comparison_metabolimics(il_dm, y_axis = "mmol/kg DM", table_title = "NMR ileum DM", 
                                            save_name = "ileum_dm", save = save)
write_tsv(il_dm_table, "tables/52_metabolomics_ileum_dm.txt")

fe_dm_table <- loop_comparison_metabolimics(fe_dm, y_axis = "mmol/kg DM", table_title = "NMR feces DM", 
                                            save_name = "feces_dm", save = save)
write_tsv(fe_dm_table, "tables/52_metabolomics_feces_dm.txt")

il_fm_table <- loop_comparison_metabolimics(il_fm, y_axis = "mmol/kg FM", table_title = "NMR ileum FM", 
                                            save_name = "ileum_fm", save = F)
write_tsv(il_fm_table, "tables/52_metabolomics_ileum_fm.txt")
fe_fm_table <- loop_comparison_metabolimics(fe_fm, y_axis = "mmol/kg FM", table_title = "NMR feces FM", 
                                            save_name = "feces_fm", save = F)
write_tsv(fe_fm_table, "tables/52_metabolomics_feces_fm.txt")