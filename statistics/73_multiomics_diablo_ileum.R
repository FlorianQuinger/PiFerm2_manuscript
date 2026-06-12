
# load image
load("temp/73_multiomics_diablo_ileum.RData")

library(here)

source(here("60_correlation_functions.R"))
source(here("70_multiomics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

#saveswitch
save = FALSE
#save = TRUE

# meta

meta <- readRDS("clean/meta1.RDS") %>%
  mutate(sampleno = str_extract(sampleid, "[0-9]{2}$")) %>%
  filter(sampleid != "103")

meta_il <- filter_ileum(meta) %>%
  arrange(sampleno) %>%
  dplyr::select(sampleno, diet, animal, period) 
meta_fa <- filter_faeces(meta) %>%
  arrange(sampleno) %>%
  dplyr::select(sampleno, diet, animal, period) %>%
  mutate(diet = as.factor(diet))

il_combined <- readRDS("temp/70_il_combined_list.RDS")

# create universe for kegg enrichment

universe <- c(colnames(il_micro$Metagenomics), colnames(il_micro$Metaproteomics)) %>%
  str_remove("[a-z]+_") 
universe <- universe[which(str_detect(universe, "K\\d{5}"))]
universe <- unique(universe)

#load data for joining


proteins_func_tax_combined <- readRDS("clean/4_function_taxa_combined_long_norm_imp_rel_abd_filtered.RDS") %>%
  inner_join(meta, by = "sampleid") %>%
  filter_ileum()

proteins_norm_imp_rel_abd_pea <- readRDS("clean/4_proteins_long_norm_imp_rel_abd.RDS") %>%
  inner_join(meta, by = "sampleid") %>%
  filter(origin %in% c("pea")) %>%
  separate(proteinid, into = c("p1", "p2", "p3"), sep = "\\|") %>%
  dplyr::select(-p1, -p2, proteinid = p3) %>%
  recalculate_rel_abd()

eggnog_mmseqs2_combined <- readRDS("clean/3_eggnog_mmseqs2_contigs_combined.RDS") %>%
  inner_join(meta, by = "sampleid")

breport_rel_abd_filtered <- readRDS("clean/3_breport_reads_long_rel_abd_filtered.RDS") %>% 
  inner_join(meta, by = "sampleid")

taxonomy_norm_imp_rel_abd <- readRDS("clean/4_taxonomy_long_norm_imp_rel_abd.RDS") %>% 
  inner_join(meta, by = "sampleid")

host_kegg_norm_imp <- readRDS("clean/4_host_kegg_long_norm_imp_rel_abd.RDS") %>%
  inner_join(meta, by = "sampleid")

# specific functions

p_plot_kegg_taxa_origin <- function(protein_taxa_df = proteins_func_tax_combined, kegg_ko_filter, taxa_level = "G") {
  filtered_df <- protein_taxa_df %>%
    filter_ileum() %>%
    filter(str_detect(kegg_ko, kegg_ko_filter)) %>%
    group_by(diet, !!sym(taxa_level), proteinid) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop")
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd, fill = !!sym(taxa_level))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors) +
    labs(title = kegg_ko_filter)
}

p_plot_taxa_abundance <- function(ktable = taxonomy_norm_imp_rel_abd, taxa_name) {
  filtered_df <- ktable %>%
    filter_ileum() %>%
    filter(rank == "G") %>%
    dplyr::select(-rank) %>%
    filter(name == taxa_name)
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey", outliers = F) +
    geom_quasirandom() +
    labs(title = taxa_name)
}

g_plot_kegg_taxa_origin <- function(gene_taxa_df = eggnog_mmseqs2_combined, kegg_ko_filter, taxa_level = "F") {
  filtered_df <- gene_taxa_df %>%
    filter_ileum() %>%
    filter(str_detect(kegg_ko, kegg_ko_filter)) %>%
    group_by(diet, sampleid, !!sym(taxa_level)) %>%
    summarise(rel_abd = sum(rel_abd, na.rm = T), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(taxa_level), values_from = "rel_abd") %>%
    pivot_longer(-c(diet, sampleid), names_to= taxa_level, values_to="rel_abd") %>%
    mutate(rel_abd = ifelse(is.na(rel_abd), 0, rel_abd)) %>%
    group_by(diet, !!sym(taxa_level)) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop")
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd, fill = !!sym(taxa_level))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors) +
    labs(title = kegg_ko_filter)
}

g_plot_taxa_abundance <- function(ktable = breport_rel_abd_filtered, taxa_name, taxa_level = "S") {
  filtered_df <- ktable %>%
    filter_ileum() %>%
    filter(rank == taxa_level) %>%
    dplyr::select(-rank) %>%
    filter(name == taxa_name) 
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey", outliers = F) +
    geom_quasirandom() +
    labs(title = taxa_name)
}

ssc_plot_kegg_abundance <- function(host_df = host_kegg_norm_imp, kegg_ko_filter) {
  filtered_df <- host_df %>%
    filter_ileum() %>%
    filter(str_detect(kegg_ko, kegg_ko_filter))
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey", outliers = F) +
    geom_quasirandom() +
    labs(title = kegg_ko_filter)
}

pea_plot_protein_abundance <- function(pea_df = proteins_norm_imp_rel_abd_pea, protein_name) {
  filtered_df <- pea_df %>%
    filter_ileum() %>%
    filter(str_detect(proteinid, protein_name))
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey", outliers = F) +
    geom_quasirandom() +
    labs(title = protein_name)
}


p_plot_kegg_for_taxa <- function(taxa_name, threshold = 1, use_names = F) {
  filtered_table <- proteins_func_tax_combined %>%
    filter(G == taxa_name) %>%
    calculate_kegg_abundance(abundance_column = "rel_abd") %>%
    group_by(sampleid) %>%
    mutate(rel_abd = rel_abd / sum(rel_abd) * 100) %>%
    ungroup() %>%
    left_join(meta, by = "sampleid") %>%
    group_by(diet, kegg_ko) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
    group_by(kegg_ko) %>%
    mutate(mean_rel_abd = mean(rel_abd)) %>%
    ungroup() %>%
    mutate(kegg_ko = ifelse(mean_rel_abd < threshold, "other", kegg_ko)) %>%
    group_by(diet, kegg_ko) %>%
    summarise(rel_abd = sum(rel_abd), .groups = "drop")
  
  if (isTRUE(use_names)) {
    filtered_table_ko <- annotate_keggs(unique(filtered_table$kegg_ko))
    
    filtered_table_kegg <- left_join(filtered_table, filtered_table_ko, by = "kegg_ko")
    
    ggplot(filtered_table_kegg, aes(x = diet, y = rel_abd, fill = reorder(name, rel_abd))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      labs(fill = "protein name", y = "relative abundance (%)", x = "")
  } else {
    ggplot(filtered_table, aes(x = diet, y = rel_abd, fill = reorder(kegg_ko, rel_abd))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      labs(fill = "KEGG ko", y = "relative abundance (%)", x = "")
  }
}

# for combined

# 1 vs 2

il_combined_1_2 <- filter_input_smccn(il_combined, meta_il, include = c("1", "2"))

cor_mean <- diablo_pairwise_correlations(il_combined_1_2[[1]])

diablo_il_combined_1_2 <- diablo_auto_model(omics_list = il_combined_1_2[[1]], meta = il_combined_1_2[[2]]$diet,
                                            design_value = cor_mean, ncomp_value = 5)

matrix_il_combined_1_2 <- plot_diablo(diablo_il_combined_1_2, cutoff_value = 0.875, save_name = "il_combined_1_2")
#matrix_il_combined_1_2 <- plot_diablo(diablo_il_combined_1_2, return_full_matrix = T)

il_combined_1_2_sub1 <- filter_submatrix(matrix_il_combined_1_2,
                                      vector = c("il_met", "il_tyr", "il_phe", "il_thr", "il_ile", "il_ala",
                                                 "il_val", "il_ser", "il_leu", "il_his", "il_asp", "il_lys",
                                                 "il_arg", "pcd_thr", "pcd_ser", "pcd_lys", "pcd_tyr", "pcd_ala",
                                                 "pcd_leu", "pcd_his", "pcd_val", "pcd_asp", "pcd_phe", "pcd_ile",
                                                 "pcd_arg", "pcd_cp", "pcd_gly",
                                                 "p_K05910", "p_K00700",
                                                 "g_K05823", 
                                                 "g_Prevotella pectinovora", "g_Prevotella sp000436035",
                                                 "g_Prevotella copri", "g_Prevotella sp900551275",
                                                 "g_Prevotella sp900554835", "g_Prevotella sp021636625",
                                                 "g_Prevotella copri_A", "g_Prevotella sp002299635",
                                                 "g_Prevotella sp004558865",
                                                 "ssc_100155945", "ssc_780428", "ssc_733607", "ssc_100514249",
                                                 "ssc_100462755", "ssc_100516959", "ssc_100157834",
                                                 "A0A9D4VFD8_PEA", "A0A9D4XAV0_PEA", "A0A9D4WNU2_PEA",
                                                 "A0A9D5APF7_PEA", "A0A9D4Y5R6_PEA", "A0A9D5AZK3_PEA",
                                                 "A0A9D4WX68_PEA", "A0A9D5BGL9_PEA", "A0A9D5AA18_PEA",
                                                 "A0A9D4YJ49_PEA", "A0A9D4VIN0_PEA", "A0A9D4W7S9_PEA",
                                                 "m_tyrosine", "m_valine", "m_phenylalanine"))

create_igraph_from_matrix(il_combined_1_2_sub1)

il_combined_1_2_sub1_sub <- filter_submatrix(matrix_il_combined_1_2,
                                             vector = c("il_met", "il_tyr", "il_phe", "il_thr", "il_ile", "il_ala",
                                                        "il_val", "il_ser", "il_leu", "il_his", "il_asp", "il_lys",
                                                        "il_arg",
                                                        "p_K05910", "p_K00700",
                                                        "g_Prevotella pectinovora", "g_Prevotella sp000436035",
                                                        "g_Prevotella copri", "g_Prevotella sp900551275",
                                                        "g_Prevotella sp900554835", "g_Prevotella sp021636625",
                                                        "g_Prevotella copri_A", "g_Prevotella sp002299635",
                                                        "g_Prevotella sp004558865",
                                                        "ssc_733607", 
                                                        "m_tyrosine", "m_valine", "m_phenylalanine"),
                                             save_name = "il_combined_1_2_sub1_sub")

create_igraph_from_matrix(il_combined_1_2_sub1_sub)

il_combined_1_2_sub1_ko <- annotate_keggs(str_remove(rownames(il_combined_1_2_sub1), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_2_sub1_ko_enriched <- enrich_keggs(filter(il_combined_1_2_sub1_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K00700")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K05910")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K10108")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K00703") # starch synthase 

p_plot_taxa_abundance(taxa_name = "Bifidobacterium")
p_plot_taxa_abundance(taxa_name = "Clostridium")
p_plot_taxa_abundance(taxa_name = "Mitsuokella")
p_plot_taxa_abundance(taxa_name = "Prevotella")

g_plot_taxa_abundance(taxa_name = "Prevotella sp000436035")
g_plot_taxa_abundance(taxa_name = "Prevotella copri_A")
g_plot_taxa_abundance(taxa_name = "Prevotella pectinovora")
g_plot_taxa_abundance(taxa_name = "Prevotella", taxa_level = "G")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K05823", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K00984", taxa_level = "P")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K02038")
#g_plot_kegg_taxa_origin(kegg_ko_filter = "K00700")
#g_plot_kegg_taxa_origin(kegg_ko_filter = "K00703")

pea_plot_protein_abundance(protein_name = "A0A9D5APF7_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4YJ49_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4VIN0_PEA")

ssc_plot_kegg_abundance(kegg_ko_filter = "780428")
ssc_plot_kegg_abundance(kegg_ko_filter = "733607")
ssc_plot_kegg_abundance(kegg_ko_filter = "100514249")
ssc_plot_kegg_abundance(kegg_ko_filter = "100157834")
ssc_plot_kegg_abundance(kegg_ko_filter = "397567")
ssc_plot_kegg_abundance(kegg_ko_filter = "100621540")
ssc_plot_kegg_abundance(kegg_ko_filter = "100462755")

p_plot_kegg_for_taxa(taxa_name = "Bifidobacterium", use_names = T, threshold = .75)
p_plot_kegg_for_taxa(taxa_name = "Mitsuokella", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Prevotella", use_names = T)


il_combined_1_2_sub1_1 <- il_combined_1_2_sub1[il_combined_1_2_sub1["ssc_780428",] != 0,
                                               il_combined_1_2_sub1[,"ssc_780428"] != 0]
create_igraph_from_matrix(il_combined_1_2_sub1_1)


il_combined_1_2_sub2 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("il_tdf",
                                                    "p_K01533", "p_K01200", "p_K01582", "p_K01585", "p_K01581",
                                                    "g_K02444", "g_K11782", "g_K02355", "g_K09949", "g_K04486",
                                                    "ssc_808504", "ssc_733685", "ssc_397113",
                                                    "A0A9D5BQA3_PEA",
                                                    "m_propionate", "m_trimethylamine"))

il_combined_1_2_sub2_ko <- annotate_keggs(str_remove(rownames(il_combined_1_2_sub2), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_2_sub2_ko_enriched <- enrich_keggs(filter(il_combined_1_2_sub2_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K01585", taxa_level = "S")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K01200")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K11782", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K09949", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K02444", taxa_level = "G")

g_plot_taxa_abundance(taxa_name = "Mitsuokella", taxa_level = "G")

ssc_plot_kegg_abundance(kegg_ko_filter = "808504")
ssc_plot_kegg_abundance(kegg_ko_filter = "397113")
ssc_plot_kegg_abundance(kegg_ko_filter = "733685")

il_combined_1_2_sub3 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("il_ip_4_1256", 
                                                    "p_K11839", "p_K09475", "p_K09476", "p_K14062", "p_K16076",
                                                    "g_K01496",
                                                    "g_Nanosynbacter sp029975625", "g_Escherichia sp004211955",
                                                    "g_Flavobacterium psychrophilum_B", "g_Escherichia sp005843885",
                                                    "g_JALHET01 sp022839525", "g_Advenella sp023423975",
                                                    "A0A9D4WQ37_PEA"))

il_combined_1_2_sub3_ko <- annotate_keggs(str_remove(rownames(il_combined_1_2_sub3), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_2_sub3_ko_enriched <- enrich_keggs(filter(il_combined_1_2_sub3_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko = "K09475")
p_plot_kegg_taxa_origin(kegg_ko = "K16076")
p_plot_kegg_taxa_origin(kegg_ko = "K09476")
p_plot_kegg_taxa_origin(kegg_ko = "K14062")

p_plot_taxa_abundance(taxa_name = "Escherichia")

g_plot_taxa_abundance(taxa_name = "Advenella sp023423975")
g_plot_taxa_abundance(taxa_name = "Escherichia sp004211955")
g_plot_taxa_abundance(taxa_name = "Escherichia sp005843885")
g_plot_taxa_abundance(taxa_name = "Nanosynbacter sp029975625")
g_plot_taxa_abundance(taxa_name = "JALHET01 sp022839525")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K01496")

g_plot_taxa_abundance(taxa_name = "Bifidobacterium", taxa_level = "G")

escherichia <- proteins_func_tax_combined %>%
  filter(G == "Escherichia") %>%
  filter(diet == 2) %>%
  group_by(proteinid, description.x, protein_name, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  mutate(sampleid = "101") %>%
  recalculate_rel_abd()

escherichia_kegg <- escherichia %>%
  calculate_kegg_abundance(abundance_column = "rel_abd")

escherichia_phosphatases <- escherichia %>%
  filter(str_detect(protein_name, "phosphatase|Phosphatase|phytase"))

escherichia_kegg_phytase <- escherichia_kegg %>%
  filter(str_detect(kegg_ko, "K01093"))

il_combined_1_2_sub4 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("m_galactose",
                                                    "p_K15582", "p_K02913", "p_K03706", "p_K02034",
                                                    "ssc_100158117"))

il_combined_1_2_sub4_ko <- annotate_keggs(str_remove(rownames(il_combined_1_2_sub4), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_2_sub4_ko_enriched <- enrich_keggs(filter(il_combined_1_2_sub4_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K02034", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K15582", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03706", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02913", taxa_level = "S")

p_plot_kegg_for_taxa(taxa_name = "Streptococcus", use_names = T)

p_plot_taxa_abundance(taxa_name = "Streptococcus")

ssc_plot_kegg_abundance(kegg_ko_filter = "100158117")

il_combined_1_2_sub5 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("pcd_ca",
                                                    "p_K00020", "p_K00042", "p_K03336", 
                                                    "g_K00012", "g_K06871",
                                                    "ssc_654405", "ssc_100152549",
                                                    "A0A9D4WM15_PEA"))

il_combined_1_2_sub5_ko <- annotate_keggs(str_remove(rownames(il_combined_1_2_sub5), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_2_sub5_ko_enriched <- enrich_keggs(filter(il_combined_1_2_sub5_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K03336")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00020")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00042")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K00012")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K06871")

ssc_plot_kegg_abundance(kegg_ko_filter = "100152549")

# 1 vs 3

il_combined_1_3 <- filter_input_smccn(il_combined, meta_il, include = c("1", "3"))

cor_mean <- diablo_pairwise_correlations(il_combined_1_3[[1]])

diablo_il_combined_1_3 <- diablo_auto_model(omics_list = il_combined_1_3[[1]], meta = il_combined_1_3[[2]]$diet,
                                            design_value = cor_mean, ncomp_value = 3)

matrix_il_combined_1_3 <- plot_diablo(diablo_il_combined_1_3, cutoff_value = 0.85, save_name = "il_combined_1_3")
diablo_il_combined_1_3$prop_expl_var
# matrix_il_combined_1_3 <- plot_diablo(diablo_il_combined_1_3, cutoff_value = 0.85, return_full_matrix = T)

g_plot_kegg_taxa_origin(kegg_ko_filter = "K02858")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K03778")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K18824")

ssc_plot_kegg_abundance(kegg_ko_filter = "780428")
ssc_plot_kegg_abundance(kegg_ko_filter = "100124377")

il_combined_1_3_sub1 <- filter_submatrix(matrix_il_combined_1_3,
                                         vector = c("enz_chy", "pcd_gly", "pcd_arg", "pcd_ala", "pcd_lys", "pcd_asp",
                                                    "pcd_his", "pcd_phe", "pcd_thr", "pcd_ile", "pcd_leu", "pcd_val",
                                                    "pcd_ser", "pcd_tyr", "enz_cara",
                                                    "il_arg", "il_lys", "il_his", "il_asp", "il_leu", "il_ala",
                                                    "il_val", "il_ile", "il_phe", "il_thr",
                                                    "p_K00186", "p_K07313", "p_K06399", "p_K21910", "p_K13923",
                                                    "g_K03778", "g_K02858", "g_K18824", "g_UBA710 sp949284615",
                                                    "ssc_100462755", "ssc_445534", "ssc_780428", "ssc_100124377",
                                                    "ssc_397604", "ssc_100526104",
                                                    "A0A9D4VJF2_PEA", "A0A9D5BG38_PEA", "A0A9D4WR70_PEA", 
                                                    "A0A9D4VKF8_PEA", "A0A9D5ANG4_PEA", "A0A9D4WDY2_PEA",
                                                    "A0A9D4YLS3_PEA", "A0A9D4X4R2_PEA", "A0A9D4VGC5_PEA",
                                                    "A0A9D4VY92_PEA", "A0A9D5AIV8_PEA", "A0A9D5AA18_PEA",
                                                    "A0A9D5GVA4_PEA", "A0A9D4YJ49_PEA", "A0A9D4Y8V2_PEA",
                                                    "A0A9D4VIN0_PEA", "A0A9D4ZZ34_PEA", "A0A9D4W1I7_PEA",
                                                    "A0A9D4XYK4_PEA", "A0A9D5AKX8_PEA", "A0A9D5AHH6_PEA", 
                                                    "m_glutamate", "m_aspartate", "m_asparagine", "m_valine",
                                                    "m_phenylalanine"))

il_combined_1_3_sub1_ko <- annotate_keggs(str_remove(rownames(il_combined_1_3_sub1), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_3_sub1_ko_enriched <- enrich_keggs(filter(il_combined_1_3_sub1_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K13923", taxa_level = "G")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K21910")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K06399", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K06398") #another sporulation protein
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01442") # bile acid metabolism
p_plot_kegg_taxa_origin(kegg_ko_filter = "K07313")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00186")

p_plot_taxa_abundance(taxa_name = "Turicibacter")
p_plot_taxa_abundance(taxa_name = "Sharpea")
p_plot_taxa_abundance(taxa_name = "Lactobacillus")
p_plot_taxa_abundance(taxa_name = "Onthomorpha")


g_plot_kegg_taxa_origin(kegg_ko_filter = "K03778", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K02858")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K18824")

g_plot_taxa_abundance(taxa_name = "UBA710 sp949284615")
g_plot_taxa_abundance(taxa_name = "Sharpea", taxa_level = "G")

ssc_plot_kegg_abundance(kegg_ko_filter = "100124377")
ssc_plot_kegg_abundance(kegg_ko_filter = "100462755")
ssc_plot_kegg_abundance(kegg_ko_filter = "397604")
ssc_plot_kegg_abundance(kegg_ko_filter = "445534")
ssc_plot_kegg_abundance(kegg_ko_filter = "780428")

pea_plot_protein_abundance(protein_name = "A0A9D5AHH6_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5BG38_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5GVA4_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4VGC5_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4ZV85_PEA")

p_plot_kegg_for_taxa(taxa_name = "Lactobacillus", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Sharpea", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Turicibacter", use_names = T)

il_combined_1_3_sub2 <- filter_submatrix(matrix_il_combined_1_3,
                                         vector = c("il_total_starch",
                                                    "p_K06305", "p_K03640", "p_K01176", "p_K02217",
                                                    "p_Basfia_A",
                                                    "g_K10542", "g_K03775", "g_K14441", "g_K03592", "g_K04085",
                                                    "ssc_397060", "ssc_100511536", "ssc_808504", "ssc_397113",
                                                    "A0A9D4VWB7_PEA", "A0A9D4VWB7_PEA", "A0A9D4ZWX1_PEA",
                                                    "A0A9D5BJS9_PEA", "A0A9D4WDI6_PEA",
                                                    "m_trimethylamine"))

il_combined_1_3_sub2_ko <- annotate_keggs(str_remove(rownames(il_combined_1_3_sub2), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_3_sub2_ko_enriched <- enrich_keggs(filter(il_combined_1_3_sub2_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K01176")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03640", taxa_level = "G")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K06305")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02217")

p_plot_taxa_abundance(taxa_name = "Basfia_A")
p_plot_taxa_abundance(taxa_name = "Turicibacter")
p_plot_taxa_abundance(taxa_name = "Escherichia")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K10542", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K03775")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K03592")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K04085")

ssc_plot_kegg_abundance(kegg_ko_filter = "100511536")
ssc_plot_kegg_abundance(kegg_ko_filter = "397113")
ssc_plot_kegg_abundance(kegg_ko_filter = "397060")
ssc_plot_kegg_abundance(kegg_ko_filter = "808504")

pea_plot_protein_abundance(protein_name = "A0A9D5BJS9_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4WDI6_PEA")

p_plot_taxa_abundance(taxa_name = "Bifidobacterium")
p_plot_taxa_abundance(taxa_name = "Pseudoscardovia")
p_plot_kegg_for_taxa(taxa_name = "Bifidobacterium", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Pseudoscardovia", use_names = T)


il_combined_1_3_sub2_1 <- il_combined_1_3_sub2[il_combined_1_3_sub2["il_total_starch",] != 0,
                                               il_combined_1_3_sub2[,"il_total_starch"] != 0]
create_igraph_from_matrix(il_combined_1_3_sub2_1)

# 1 vs 4

il_combined_1_4 <- filter_input_smccn(il_combined, meta_il, include = c("1", "4"))

cor_mean <- diablo_pairwise_correlations(il_combined_1_4[[1]])

diablo_il_combined_1_4 <- diablo_auto_model(omics_list = il_combined_1_4[[1]], meta = il_combined_1_4[[2]]$diet,
                                        design_value = cor_mean, ncomp_value = 4)

matrix_il_combined_1_4 <- plot_diablo(diablo_il_combined_1_4, cutoff_value = 0.8, save_name = "il_combined_1_4")
diablo_il_combined_1_4$prop_expl_var
#matrix_il_combined_1_4 <- plot_diablo(diablo_il_combined_1_4, return_full_matrix = T)

g_plot_kegg_taxa_origin(kegg_ko_filter = "K12984")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K05934")

p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03074")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03072")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K12257", taxa_level = "G")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K07699", taxa_level = "S")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03413")

ssc_plot_kegg_abundance(kegg_ko_filter = "100620246")
ssc_plot_kegg_abundance(kegg_ko_filter = "396919")

il_combined_1_4_sub1 <- filter_submatrix(matrix_il_combined_1_4,
                                      vector = c("pcd_phe", "pcd_pro", "il_pro", "pcd_gly", "il_gly",
                                                 "p_K03074", "p_K03072", "p_K03413", "p_K07699", "p_K12257",
                                                 "g_K02169", "g_K05934", "g_K12984", 
                                                 "g_Bifidobacterium boum",
                                                 "ssc_396919", "ssc_100620246", "ssc_396807", "ssc_396766", 
                                                 "ssc_397197", "ssc_445461", "ssc_100514249", "ssc_100157995",
                                                 "ssc_100158244", "ssc_397316", 
                                                 "A0A9D4ZV85_PEA", "A0A9D4VX07_PEA", "A0A9D4XPB1_PEA",
                                                 "A0A9D5A4B1_PEA", "A0A9D5A4J7_PEA",
                                                 "m_glutamine", "m_phenylalanine", "m_valine", "m_aspartate", 
                                                 "m_asparagine", "m_aspartate", "m_glutamate", "m_trimethylamine"))

il_combined_1_4_sub1_ko <- annotate_keggs(str_remove(rownames(il_combined_1_4_sub1), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_4_sub1_ko_enriched <- enrich_keggs(filter(il_combined_1_4_sub1_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03413")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K07699")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03072")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03074", taxa_level = "P")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K12257")

p_plot_taxa_abundance(taxa_name = "Ruminococcus")
p_plot_taxa_abundance(taxa_name = "HUN007")
p_plot_taxa_abundance(taxa_name = "Onthomorpha")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K02169")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K12984", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K05934")

g_plot_taxa_abundance(taxa_name = "Bifidobacterium boum")

ssc_plot_kegg_abundance(kegg_ko_filter = "397197")
ssc_plot_kegg_abundance(kegg_ko_filter = "396766")
ssc_plot_kegg_abundance(kegg_ko_filter = "396919")
ssc_plot_kegg_abundance(kegg_ko_filter = "396807")
ssc_plot_kegg_abundance(kegg_ko_filter = "445461")
ssc_plot_kegg_abundance(kegg_ko_filter = "100158244")
ssc_plot_kegg_abundance(kegg_ko_filter = "397316")
ssc_plot_kegg_abundance(kegg_ko_filter = "397316")

pea_plot_protein_abundance(protein_name = "A0A9D4XPB1_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5A4B1_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5A4J7_PEA")

p_plot_kegg_for_taxa(taxa_name = "Prevotella", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Onthomorpha", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Mitsuokella", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Limisoma", use_names = T)

il_combined_1_4_sub2 <- filter_submatrix(matrix_il_combined_1_4,
                                         vector = c("il_total_starch", "il_insp6", "il_ca", "pcd_ca", "pcd_p",
                                                    "p_K12686", "p_K03281", "p_K07571", "p_K02372", "p_Enterocloster",
                                                    "g_K00863", "g_K01961", "g_K03718",
                                                    "g_Hominenteromicrobium mulieris", "g_Bacteroides xylanisolvens",
                                                    "g_Bariatricus sp004560705", 
                                                    "ssc_100521982", "ssc_100624628", "ssc_397602", "ssc_396921",
                                                    "A0A9D4XTF6_PEA", "A0A9D5BQA3_PEA", "A0A9D4XKQ8_PEA",
                                                    "A0A9D5B8M5_PEA",
                                                    "m_trimethylamine", "m_propionate"))

il_combined_1_4_sub2_ko <- annotate_keggs(str_remove(rownames(il_combined_1_4_sub2), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_4_sub2_ko_enriched <- enrich_keggs(filter(il_combined_1_4_sub2_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K02372")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01686")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K12686")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K07571")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03281")

p_plot_taxa_abundance(taxa_name = "Enterocloster")
p_plot_taxa_abundance(taxa_name = "Mitsuokella")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K00863")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K01961")

g_plot_taxa_abundance(taxa_name = "Bacteroides xylanisolvens")
g_plot_taxa_abundance(taxa_name = "Hominenteromicrobium mulieris")

ssc_plot_kegg_abundance(kegg_ko_filter = "397602")
ssc_plot_kegg_abundance(kegg_ko_filter = "100624628")
ssc_plot_kegg_abundance(kegg_ko_filter = "396921")
ssc_plot_kegg_abundance(kegg_ko_filter = "100521982")


il_combined_1_4_sub3 <- filter_submatrix(matrix_il_combined_1_4,
                                         vector = c("enz_try", "enz_chy", "enz_amy", "il_p", "il_total_starch", 
                                                    "p_K17810", "p_K01493", "p_K05910", "p_K03218",
                                                    "g_K03980", "g_K06956", "g_K02502", "g_K03718",
                                                    "g_Fusobacterium gastrosuis", 
                                                    "ssc_100155038", "ssc_808504", "ssc_397113",
                                                    "A0A9D5BN09_PEA", "A0A9D4YP35_PEA", "A0A9D4WQC0_PEA",
                                                    "A0A9D4VWB7_PEA", "A0A9D4YG61_PEA",
                                                    "m_lactate", "m_tyrosine"))

il_combined_1_4_sub3_ko <- annotate_keggs(str_remove(rownames(il_combined_1_4_sub3), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_4_sub3_ko_enriched <- enrich_keggs(filter(il_combined_1_4_sub3_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K01493")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K05910", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17810", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03218", taxa_level = "S")

ssc_plot_kegg_abundance(kegg_ko_filter = "397113")

ssc_plot_kegg_abundance(kegg_ko_filter = "100155038")

g_plot_taxa_abundance(taxa_name = "Fusobacterium gastrosuis")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K03980")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K02502")

pea_plot_protein_abundance(protein_name = "A0A9D5BN09_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4YP35_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4WQC0_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4VWB7_PEA")


il_combined_1_4_sub3_1 <- il_combined_1_4_sub3[il_combined_1_4_sub3["ssc_397113",] != 0,
                                               il_combined_1_4_sub3[,"ssc_397113"] != 0]
create_igraph_from_matrix(il_combined_1_4_sub3_1)


il_combined_1_4_sub4 <- filter_submatrix(matrix_il_combined_1_4,
                                         vector = c("p_K00854", "p_K01738", "p_K19168", "p_K00962",
                                                    "g_K00873", "g_K14261", "g_K06074",
                                                    "ssc_595116", 
                                                    "A0A9D4VQB0_PEA",
                                                    "m_butyrate"))

il_combined_1_4_sub4_ko <- annotate_keggs(str_remove(rownames(il_combined_1_4_sub4), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_4_sub4_ko_enriched <- enrich_keggs(filter(il_combined_1_4_sub4_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00854")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K19168")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00962")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01738")

ssc_plot_kegg_abundance(kegg_ko_filter = "595116")

pea_plot_protein_abundance(protein_name = "A0A9D4VQB0_PEA")

il_combined_1_4_sub5 <- filter_submatrix(matrix_il_combined_1_4,
                                         vector = c("p_K06942",
                                                    "ssc_100526104",
                                                    "A0A9D5B868_PEA", "A0A9D4W9R0_PEA"))

il_combined_1_4_sub5_ko <- annotate_keggs(str_remove(rownames(il_combined_1_4_sub5), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
il_combined_1_4_sub5_ko_enriched <- enrich_keggs(filter(il_combined_1_4_sub5_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

ssc_plot_kegg_abundance(kegg_ko_filter = "100526104")

pea_plot_protein_abundance(protein_name = "A0A9D4W9R0_PEA")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K06942")

# save environment
save.image("temp/73_multiomics_diablo_ileum.RData")

