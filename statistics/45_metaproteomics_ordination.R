library(here)
library(vegan)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

set.seed(1112)

#save switch
#save = TRUE
save = FALSE

# load meta

meta <- readRDS("clean/meta1.RDS")  %>% filter(sampleid != "103")

# load proteins file with rel_abd

proteins_raw <- readRDS("clean/4_proteins_long_raw_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
proteins_norm_imp <- readRDS("clean/4_proteins_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

taxa_raw <- readRDS("clean/4_taxa_long_raw_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
taxa_norm_imp <- readRDS("clean/4_taxa_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

kegg_raw <- readRDS("clean/4_kegg_long_raw_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
kegg_norm_imp <- readRDS("clean/4_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

host_kegg_raw <- readRDS("clean/4_host_kegg_long_raw_rel_abd_filtered.RDS") %>% filter(sampleid != "103")
host_kegg_norm_imp <- readRDS("clean/4_host_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")


# pcoa for proteins

do_ordination(proteins_raw, title = "Proteingroups raw", save_name = "proteingroups_raw", save = save)
do_ordination(proteins_norm_imp, title = "Proteingroups imputed", save_name = "proteingroups_imp", save = save)
do_ordination(proteins_raw, region = "ileum") # sig
do_ordination(proteins_norm_imp, region = "ileum", title = "Proteingroups ileum imputed", save_name = "proteingroups_ileum_imp", save = save) #sig
do_ordination(proteins_raw, region = "faeces") #sig
do_ordination(proteins_norm_imp, region = "faeces", title = "Proteingroups faeces imputed", save_name = "proteingroups_faeces_imp", save = save) #sig

# filter for host proteins

proteins_raw_host <- proteins_raw %>%
  filter(origin == "pig") %>%
  recalculate_rel_abd()
proteins_norm_imp_host <- proteins_norm_imp %>%
  filter(origin == "pig") %>%
  recalculate_rel_abd()

do_ordination(proteins_raw_host)
do_ordination(proteins_norm_imp_host)
do_ordination(proteins_raw_host, region = "ileum", second_indicator = "animal") #sig
do_ordination(proteins_norm_imp_host, region = "ileum", title = "Proteingroups host ileum", save_name = "proteingroups_host_ileum_imp", save = T) #sig
do_ordination(proteins_raw_host, region = "faeces", second_indicator = "animal") #sig
do_ordination(proteins_norm_imp_host, region = "faeces", title = "Proteingroups host faeces", save_name = "proteingroups_host_faeces_imp", save = T) #sig

# filter for bacterial proteins

proteins_raw_micro <- proteins_raw %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()
proteins_norm_imp_micro <- proteins_norm_imp %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()

do_ordination(proteins_raw_micro)
do_ordination(proteins_norm_imp_micro)
do_ordination(proteins_raw_micro, region = "ileum")
do_ordination(proteins_norm_imp_micro, region = "ileum", title = "Proteingroups micro ileum", save_name = "proteingroups_micro_ileum_imp", save = T) 
do_ordination(proteins_raw_micro, region = "faeces", second_indicator = "animal") #sig but not pairwise
do_ordination(proteins_norm_imp_micro, region = "faeces", title = "Proteingroups micro faeces", save_name = "proteingroups_micro_faeces_imp", save = T) #sig but not pairwise

# for taxa

do_ordination(taxa_raw)
do_ordination(taxa_norm_imp)
do_ordination(taxa_raw, region = "ileum")
do_ordination(taxa_norm_imp, region = "ileum", title = "Taxa ileum imputed", save_name = "taxa_ileum_imp", save = T)
do_ordination(taxa_raw, region = "faeces") #sig but not pairwise
do_ordination(taxa_norm_imp, region = "faeces", title = "Taxa faeces imputed", save_name = "taxa_faeces_imp", save = T) #sig but not pairwise

# for kegg functions

do_ordination(kegg_raw)
do_ordination(kegg_norm_imp)
do_ordination(kegg_raw, region = "ileum")
do_ordination(kegg_norm_imp, region = "ileum", title = "Kegg micro ileum", save_name = "kegg_micro_ileum_imp", save = T)
do_ordination(kegg_raw, region = "faeces") #sig but not pairwise
do_ordination(kegg_norm_imp, region = "faeces", title = "Kegg micro faeces", save_name = "kegg_micro_faeces_imp", save = T) #sig  but not pairwise

# for host kegg functions

do_ordination(host_kegg_raw)
do_ordination(host_kegg_norm_imp)
do_ordination(host_kegg_raw, region = "ileum") 
do_ordination(host_kegg_norm_imp, region = "ileum", title = "Kegg host ileum", save_name = "kegg_host_ileum_imp", save = T) 
do_ordination(host_kegg_raw, region = "faeces") #sig 
do_ordination(host_kegg_norm_imp, region = "faeces", title = "Kegg host faeces", save_name = "kegg_host_faeces_imp", save = T) #sig

# for feed proteins

proteins_raw_feed <- proteins_raw %>%
  filter(origin %in% c("soy", "pea", "rape", "barley", "wheat")) %>%
  recalculate_rel_abd()
proteins_norm_imp_feed <- proteins_norm_imp %>%
  filter(origin %in% c("soy", "pea", "rape", "barley", "wheat")) %>%
  recalculate_rel_abd()

do_ordination(proteins_raw_feed)
do_ordination(proteins_norm_imp_feed)
do_ordination(proteins_raw_feed, region = "ileum") # sig
do_ordination(proteins_norm_imp_feed, region = "ileum", title = "Proteingroups feed ileum", save_name = "proteingroups_feed_ileum_imp", save = T) # sig
do_ordination(proteins_raw_feed, region = "faeces") #sig
do_ordination(proteins_norm_imp_feed, region = "faeces", title = "Proteingroups feed faeces", save_name = "proteingroups_feed_faeces_imp", save = T) #sig 

# for only pea proteins

proteins_norm_imp_pea <- proteins_norm_imp %>%
  filter(origin %in% c("pea")) %>%
  recalculate_rel_abd()

do_ordination(proteins_norm_imp_pea, region = "ileum")
do_ordination(proteins_norm_imp_pea, region = "faeces")
