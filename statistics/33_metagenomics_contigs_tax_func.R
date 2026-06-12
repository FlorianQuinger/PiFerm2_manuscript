library(here)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

#save switch
#save = TRUE
save = FALSE

# load meta

meta <- readRDS("clean/meta1.RDS") %>% filter(sampleid != "103")

# load function files

#eggnog_reads <- readRDS("clean/3_eggnog_contigs_long_reads.RDS")
eggnog_reads_filtered <- readRDS("clean/3_eggnog_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")

#eggnog_rel_abd <- readRDS("clean/3_eggnog_contigs_long_rel_abd.RDS")
eggnog_rel_abd_filtered <- readRDS("clean/3_eggnog_contigs_long_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

eggnog <- readRDS("clean/3_eggnog_contigs_long_raw.RDS") %>%
  dplyr::select(-c(sampleid, rel_abd, reads, gene, evalue, score)) %>%
  distinct()

#kegg_reads <- readRDS("clean/3_kegg_contigs_long_reads.RDS")
kegg_reads_filtered <- readRDS("clean/3_kegg_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")

#kegg_rel_abd <- readRDS("clean/3_kegg_contigs_long_rel_abd.RDS")
kegg_rel_abd_filtered <- readRDS("clean/3_kegg_contigs_long_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

#go_reads <- readRDS("clean/3_go_contigs_long_reads.RDS")
go_reads_filtered <- readRDS("clean/3_go_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")

#go_rel_abd <- readRDS("clean/3_go_contigs_long_rel_abd.RDS")
go_rel_abd_filtered <- readRDS("clean/3_go_contigs_long_rel_abd_filtered.RDS") %>% filter(sampleid != "103")


cog_reads_filtered <- readRDS("clean/3_cog_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")
cog_rel_abd_filtered <- readRDS("clean/3_cog_contigs_long_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

cazy_reads_filtered <- readRDS("clean/3_cazy_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")
cazy_rel_abd_filtered <- readRDS("clean/3_cazy_contigs_long_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

# Ordination based on functions

# for all functions

do_ordination(eggnog_rel_abd_filtered, region = "ileum", title = "Functions contigs ileum", 
              save_name = "func_ileum_contigs", save = save)

do_ordination(eggnog_rel_abd_filtered, region = "faeces", title = "Functions contigs faeces", 
              save_name = "func_faeces_contigs", save = save)

# for keggs

do_ordination(kegg_rel_abd_filtered, region = "ileum", title = "KEGG contigs ileum", 
              save_name = "kegg_ileum_contigs", save = save)

do_ordination(kegg_rel_abd_filtered, region = "faeces", title = "KEGG contigs faeces", 
              save_name = "kegg_faeces_contigs", save = save)

# for go's

do_ordination(go_rel_abd_filtered, region = "ileum", title = "GO contigs ileum", 
              save_name = "go_ileum_contigs", save = save)

do_ordination(go_rel_abd_filtered, region = "faeces", title = "GO contigs faeces", 
              save_name = "go_faeces_contigs", save = save)

# for cog's

do_ordination(cog_rel_abd_filtered, region = "ileum", title = "COG contigs ileum", 
              save_name = "cog_ileum_contigs", save = save)

do_ordination(cog_rel_abd_filtered, region = "faeces", title = "COG contigs faeces", 
              save_name = "cog_faeces_contigs", save = save)

# for cog's

do_ordination(cazy_rel_abd_filtered, region = "ileum", title = "CAZY contigs ileum", 
              save_name = "cazy_ileum_contigs", save = save)

do_ordination(cazy_rel_abd_filtered, region = "faeces", title = "CAZY contigs faeces", 
              save_name = "cazy_faeces_contigs", save = save)

# differential abundance of functions

il_eggnog <- perform_edger_and_plot(eggnog_reads_filtered, selected_matrix = "ileum", save = save,
                                    abundance_column = "reads", save_name = "eggnog_contigs")

fa_eggnog <- perform_edger_and_plot(eggnog_reads_filtered,selected_matrix = "faeces", save = save,
                                    abundance_column = "reads", save_name = "eggnog_contigs")

il_eggnog_func <- add_functions(results_object = il_eggnog, functions_object = eggnog, join_by = "seed_ortholog")

fa_eggnog_func <- add_functions(results_object = fa_eggnog, functions_object = eggnog, join_by = "seed_ortholog")

il_eggnog_kegg <- transform_to_function(il_eggnog_func, function_name = "kegg_ko")

il_eggnog_kegg <- transform_to_function(fa_eggnog_func, function_name = "kegg_ko")

# enrichment analysis with translated functions

enrich_kegg(il_eggnog_kegg, save_name = "33_enrichment_eggnog_kegg_contigs_ileum", save = save)

enrich_kegg(fa_eggnog_kegg, save_name = "33_enrichment_eggnog_kegg_contigs_faeces", save = save)

# differential abundance of keggs

il_kegg <- perform_edger_and_plot(kegg_reads_filtered, selected_matrix = "ileum", save = save,
                                  abundance_column = "reads", save_name = "kegg_contigs")

enrich_kegg(il_kegg, save_name = "33_enrichment_kegg_contigs_ileum", save = save)

fa_kegg <- perform_edger_and_plot(kegg_reads_filtered, selected_matrix = "faeces", save = save,
                                  abundance_column = "reads", save_name = "kegg_contigs")

enrich_kegg(fa_kegg, save_name = "33_enrichment_kegg_contigs_faeces", save = save)

source("70_multiomics_functions.R")

add_kegg_annotation <- function(edger_object) {
  for (i in 1:length(edger_object)) {
    table_sig <- edger_object[[i]][["table_sig"]] 
    annotated_keggs <- table_sig %>% 
      filter(sig == "sig") %>%
      pull(kegg_ko) %>%
      annotate_keggs()
    
    table_sig <- table_sig %>%
      left_join(annotated_keggs, by = "kegg_ko")
    
    edger_object[[i]][["table_sig"]] <- table_sig
  }
  return(edger_object)
}

il_kegg_annotated <- add_kegg_annotation(edger_object = il_kegg)

fa_kegg_annotated <- add_kegg_annotation(edger_object = fa_kegg)

# differential abundance of  gos

il_go <- perform_edger_and_plot(go_reads_filtered, selected_matrix = "ileum", save = save,
                                  abundance_column = "reads", save_name = "go_contigs")

fa_go <- perform_edger_and_plot(go_reads_filtered, selected_matrix = "faeces", save = save,
                                  abundance_column = "reads", save_name = "go_contigs")

# differential abundance of  cogs

il_cog <- perform_edger_and_plot(cog_reads_filtered, selected_matrix = "ileum", save = save,
                                abundance_column = "reads", save_name = "cog_contigs")

fa_cog <- perform_edger_and_plot(cog_reads_filtered, selected_matrix = "faeces", save = save,
                                abundance_column = "reads", save_name = "cog_contigs")

# enrichment analysis of cogs

enrich_cog(il_cog, save_name = "33_enrichment_cog_ileum") 
enrich_cog(fa_cog, save_name = "33_enrichment_cog_faeces")

# differential abundance of  cazys

il_cazy <- perform_edger_and_plot(cazy_reads_filtered, selected_matrix = "ileum", save = save,
                                 abundance_column = "reads", save_name = "cazy_contigs") # no significant differences

fa_cazy <- perform_edger_and_plot(cazy_reads_filtered, selected_matrix = "faeces", save = save,
                                 abundance_column = "reads", save_name = "cazy_contigs")

# extract significant functions

kegg_il_sig <- create_venn_for_sig(il_kegg, kegg_rel_abd_filtered, 
                                   selected_matrix = "ileal digesta", save_name = "kegg", save = save)
kegg_fa_sig <- create_venn_for_sig(fa_kegg, kegg_rel_abd_filtered, 
                                   selected_matrix = "faeces", save_name = "kegg", save = save)

go_il_sig <- create_venn_for_sig(il_go, go_rel_abd_filtered,
                                 selected_matrix = "ileal digesta", save_name = "go", save = save)
go_fa_sig <- create_venn_for_sig(fa_go, go_rel_abd_filtered,
                                 selected_matrix = "faeces", save_name = "go", save = save) 

cog_il_sig <- create_venn_for_sig(il_cog, cog_rel_abd_filtered,
                                  selected_matrix = "ileal digesta", save_name = "cog", save = save)
cog_fa_sig <- create_venn_for_sig(fa_cog, cog_rel_abd_filtered,
                                  selected_matrix = "faeces", save_name = "cog", save = save)

#cazy_il_sig <- create_venn_for_sig(il_cazy, cazy_rel_abd_filtered,
#                                   selected_matrix = "ileal digesta", save_name = "cazy", save = save)
cazy_fa_sig <- create_venn_for_sig(fa_cazy, cazy_rel_abd_filtered,
                                   selected_matrix = "faeces", save_name = "cazy", save = save)

# create heatmaps for functions

create_heatmap_from_edger(il_kegg, kegg_il_sig, heatmap_y = "kegg_ko", save = save,
                          save_name = "33_heatmap_kegg_ileum")
create_heatmap_from_edger(fa_kegg, kegg_fa_sig, heatmap_y = "kegg_ko", save = save,
                          save_name = "33_heatmap_kegg_faeces")

create_heatmap_from_edger(il_go, go_il_sig, heatmap_y = "go", save = save,
                          save_name = "33_heatmap_go_ileum")
create_heatmap_from_edger(fa_go, go_fa_sig, heatmap_y = "go", save = save,
                          save_name = "33_heatmap_go_faeces")

create_heatmap_from_edger(il_cog, cog_il_sig, heatmap_y = "cog_accession", save = save,
                          save_name = "33_heatmap_cog_ileum")
create_heatmap_from_edger(fa_cog, cog_fa_sig, heatmap_y = "cog_accession", save = save,
                          save_name = "33_heatmap_cog_faeces")

create_heatmap_from_edger(fa_cazy, cazy_fa_sig, heatmap_y = "cazy", save = save,
                          save_name = "33_heatmap_cazy_faeces")

