
# load image
load("temp/74_multiomics_diablo_faeces.RData")

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

fa_combined <- readRDS("temp/70_fa_combined_list.RDS")

# create universe for kegg enrichment

universe <- c(colnames(fa_micro$Metagenomics), colnames(fa_micro$Metaproteomics)) %>%
  str_remove("[a-z]+_") 
universe <- universe[which(str_detect(universe, "K\\d{5}"))]
universe <- unique(universe)

#load data for joining

proteins_func_tax_combined <- readRDS("clean/4_function_taxa_combined_long_norm_imp_rel_abd_filtered.RDS") %>%
  inner_join(meta, by = "sampleid") %>%
  filter_faeces()

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
    filter_faeces() %>%
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
    filter_faeces() %>%
    filter(rank == "G") %>%
    dplyr::select(-rank) %>%
    filter(name == taxa_name)
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey") +
    geom_quasirandom() +
    labs(title = taxa_name)
}

g_plot_kegg_taxa_origin <- function(gene_taxa_df = eggnog_mmseqs2_combined, kegg_ko_filter, taxa_level = "F") {
  filtered_df <- gene_taxa_df %>%
    filter_faeces() %>%
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
    filter_faeces() %>%
    filter(rank == taxa_level) %>%
    dplyr::select(-rank) %>%
    filter(name == taxa_name) 
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey") +
    geom_quasirandom() +
    labs(title = taxa_name)
}

ssc_plot_kegg_abundance <- function(host_df = host_kegg_norm_imp, kegg_ko_filter) {
  filtered_df <- host_df %>%
    filter_faeces() %>%
    filter(str_detect(kegg_ko, kegg_ko_filter))
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey") +
    geom_quasirandom() +
    labs(title = kegg_ko_filter)
}

pea_plot_protein_abundance <- function(pea_df = proteins_norm_imp_rel_abd_pea, protein_name) {
  filtered_df <- pea_df %>%
    filter_faeces() %>%
    filter(str_detect(proteinid, protein_name))
  
  ggplot(filtered_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(fill = "grey") +
    geom_quasirandom() +
    labs(title = protein_name)
}

p_plot_kegg_for_taxa <- function(taxa_name, threshold = 1, use_names = F, taxa_level = "G") {
  filtered_table <- proteins_func_tax_combined %>%
    filter(!!sym(taxa_level) == taxa_name) %>%
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
    mutate(kegg_ko = ifelse(mean_rel_abd <= threshold, "other", kegg_ko)) %>%
    group_by(diet, kegg_ko) %>%
    summarise(rel_abd = sum(rel_abd), .groups = "drop")
  
  if (isTRUE(use_names)) {
    filtered_table_ko <- annotate_keggs(unique(filtered_table$kegg_ko))
    
    filtered_table_kegg <- left_join(filtered_table, filtered_table_ko, by = "kegg_ko")
    
    p <- ggplot(filtered_table_kegg, aes(x = diet, y = rel_abd, fill = reorder(name, rel_abd))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      labs(fill = "protein name", y = "relative abundance (%)", x = "")
    print(p)
    return(filtered_table_kegg)
  } else {
    p <- ggplot(filtered_table, aes(x = diet, y = rel_abd, fill = reorder(kegg_ko, rel_abd))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      labs(fill = "KEGG ko", y = "relative abundance (%)", x = "")
    print(p)
  }
}

# diablo for combined data 

# 1 vs 2

fa_combined_1_2 <- filter_input_smccn(fa_combined, meta_fa, include = c("1", "2"))

cor_mean <- diablo_pairwise_correlations(fa_combined_1_2[[1]])

diablo_fa_combined_1_2 <- diablo_auto_model(omics_list = fa_combined_1_2[[1]], meta = fa_combined_1_2[[2]]$diet,
                                            design_value = cor_mean, ncomp_value = 3)

matrix_fa_combined_1_2 <- plot_diablo(diablo_fa_combined_1_2, cutoff = 0.75, save_name = "fa_combined_1_2")
diablo_fa_combined_1_2$prop_expl_var
#matrix_fa_combined_1_2 <- plot_diablo(diablo_fa_combined_1_2, return_full_matrix = T)
#network_fa_combined_1_2 <- plot_diablo(diablo_fa_combined_1_2, cutoff = 0.5, return_network = T)

g_plot_kegg_taxa_origin(kegg_ko_filter = "K01267", taxa_level = "F")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K02535", taxa_level = "F")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K16363", taxa_level = "F")

p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K02243")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K02652")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03544")

p_plot_taxa_abundance(taxa_name = "Phascolarctobacterium")

ssc_plot_kegg_abundance(kegg_ko_filter = "396685")
ssc_plot_kegg_abundance(kegg_ko_filter = "100521982")
ssc_plot_kegg_abundance(kegg_ko_filter = "100312960")
ssc_plot_kegg_abundance(kegg_ko_filter = "406192")

fa_combined_1_2_sub1 <- filter_submatrix(matrix_fa_combined_1_2,
                                         vector = c("hg_ge", "hg_cp", "fa_p",
                                                    "p_K00440", "p_K02032", "p_K07386", "p_K01622", "p_K04088",
                                                    "p_K04087", "p_K10192", "p_K03837",
                                                    "p_RGIG1693", "p_Hominicoprocola",
                                                    "g_K07068", "g_K07701", "g_K17319", "g_K01676",
                                                    "g_CAG-103 sp905215475",
                                                    "ssc_100627480", "ssc_397140", "ssc_397419",
                                                    "A0A9D5B2R9_PEA", "A0A9D4VRF6_PEA", "A0A068LJH6_PEA",
                                                    "m_3-phenylpropionate", "m_valerate", "m_propionate"))

fa_combined_1_2_sub1_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_2_sub1), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_2_sub1_ko_enriched <- enrich_keggs(filter(fa_combined_1_2_sub1_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K04088")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K01622")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K07386")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K02032")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K03837")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K10192")
p_plot_kegg_taxa_origin(proteins_func_tax_combined, kegg_ko_filter = "K00440")

p_plot_taxa_abundance(taxa_name = "Hominicoprocola")
p_plot_taxa_abundance(taxa_name = "Cryptobacteroides")
p_plot_taxa_abundance(taxa_name = "Faecousia")
p_plot_taxa_abundance(taxa_name = "Methanobrevibacter_A")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K17319")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K01676")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K07068")

ssc_plot_kegg_abundance(kegg_ko_filter = "397419")
ssc_plot_kegg_abundance(kegg_ko_filter = "397140")

g_plot_taxa_abundance(taxa_name = "CAG-103 sp905215475")

p_plot_kegg_for_taxa(taxa_name = "Methanobrevibacter_A", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Hominicoprocola", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Faecousia", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Treponema_D", use_names = T)

hominicoprocola <- proteins_func_tax_combined %>%
  filter(G == "Hominicoprocola") %>%
  filter(diet == 2) %>%
  group_by(proteinid, description.x, protein_name, kegg_ko) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
  mutate(sampleid = "101") %>%
  recalculate_rel_abd()

fa_combined_1_2_sub2 <- filter_submatrix(matrix_fa_combined_1_2,
                                         vector = c("fa_ti", "fa_cp",
                                                    "p_K07335", "p_K07337", "p_K03832", "p_K00348", 
                                                    "p_PeH17", "p_UBA4363",
                                                    "g_K01846", "g_K01442", "g_K07016",
                                                    "ssc_397192", "ssc_100522855", "ssc_100154047",
                                                    "A0A9D4W4Y0_PEA", "A0A9D4WKA7_PEA", "D3VND9_PEA", "Q9M3X6_PEA",
                                                    "A0A9D4XWP6_PEA",
                                                    "m_valine", "m_phenylalanine", "m_methionine", "m_tyrosine"))
create_igraph_from_matrix(fa_combined_1_2_sub2)

fa_combined_1_2_sub2_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_2_sub2), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_2_sub2_ko_enriched <- enrich_keggs(filter(fa_combined_1_2_sub2_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K07337", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K07335", taxa_level = "G" )
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03832")

p_plot_taxa_abundance(taxa_name = "Treponema_D")
p_plot_taxa_abundance(taxa_name = "PeH17")
p_plot_taxa_abundance(taxa_name = "Prevotella")
p_plot_taxa_abundance(taxa_name = "UBA4363")

ssc_plot_kegg_abundance(kegg_ko_filter = "100522855")
ssc_plot_kegg_abundance(kegg_ko_filter = "100154047")
ssc_plot_kegg_abundance(kegg_ko_filter = "397192")

pea_plot_protein_abundance(protein_name = "A0A9D4XWP6_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4W4Y0_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4WKA7_PEA")
pea_plot_protein_abundance(protein_name = "D3VND9_PEA")
pea_plot_protein_abundance(protein_name = "Q9M3X6_PEA")

p_plot_kegg_for_taxa(taxa_name = "PeH17", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Treponema_D", use_names = T)

fa_combined_1_2_sub2_1 <- fa_combined_1_2_sub2[fa_combined_1_2_sub2["p_PeH17",] != 0,
                                               fa_combined_1_2_sub2[,"p_PeH17"] != 0]
create_igraph_from_matrix(fa_combined_1_2_sub2_1)


fa_combined_1_2_sub3 <- filter_submatrix(matrix_fa_combined_1_2,
                                         vector = c("enz_try", "fa_insp6", "hg_insp6",
                                                    "p_K02652", "p_K02243", "p_K03544", "p_K22339",
                                                    "p_Phascolarctobacterium", 
                                                    "g_K01267", "g_K02535", "g_K11732", "g_K16363", 
                                                    "g_Ruminococcus_C sp937915125",
                                                    "ssc_100521982", "ssc_406192", "ssc_733683",
                                                    "ssc_396685", "ssc_100312960",
                                                    "A0A9D4VYD1_PEA", "A0A9D5BH02_PEA", "A0A9D4Y6S0_PEA",
                                                    "A0A9D4VY92_PEA", "A0A9D4X2N1_PEA",
                                                    "m_succinate", "m_acetate"))

fa_combined_1_2_sub3_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_2_sub3), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_2_sub3_ko_enriched <- enrich_keggs(filter(fa_combined_1_2_sub3_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K02243")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02652")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03544")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K22339", taxa_level = "S")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K01938") # other formate tetrahydrofolate ligase from faecousia
#hydrogenases
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00436")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18005")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18006")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18007")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17992")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17993")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17994")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18330")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18331")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18332")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K17997")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17998")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17999")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17993")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17994")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00437")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18008")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K05922")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K05927")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00532")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00533")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00534")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K06441")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18016")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18017")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18023")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00440")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00441")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00442")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K13942")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K14068")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K14069")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K14070")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17995")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K17996")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K06281")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K06282")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K23548")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K23549")
# other pathways
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00399")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00239")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K16950")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00370")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K03385")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K14138")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00198")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02117")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00425")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00937") # polyp kinase

p_plot_taxa_abundance(taxa_name = "Faecousia")
p_plot_taxa_abundance(taxa_name = "Phascolarctobacterium")
p_plot_taxa_abundance(taxa_name = "Ruminococcus")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K01267")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K16363")

g_plot_taxa_abundance(taxa_name = "Ruminococcus_C sp937915125")
g_plot_taxa_abundance(taxa_name = "Megasphaeraceae", taxa_level = "F")

ssc_plot_kegg_abundance(kegg_ko_filter = "100521982")
ssc_plot_kegg_abundance(kegg_ko_filter = "396685")

p_plot_kegg_for_taxa(taxa_name = "Ruminococcus", use_names = T, threshold = 0.75)
p_plot_kegg_for_taxa(taxa_name = "Faecousia", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Phascolarctobacterium", use_names = T)
p <- p_plot_kegg_for_taxa(taxa_name = "Megasphaera", use_names = T, threshold = 0)

# 1 vs 3

fa_combined_1_3 <- filter_input_smccn(fa_combined, meta_fa, include = c("1", "3"))

cor_mean <- diablo_pairwise_correlations(fa_combined_1_3[[1]])

diablo_fa_combined_1_3 <- diablo_auto_model(omics_list = fa_combined_1_3[[1]], meta = fa_combined_1_3[[2]]$diet,
                                            design_value = cor_mean, ncomp_value = 2)

matrix_fa_combined_1_3 <- plot_diablo(diablo_fa_combined_1_3, cutoff_value = .8, save_name = "fa_combined_1_3")
#matrix_fa_combined_1_3 <- plot_diablo(diablo_fa_combined_1_3, return_full_matrix = T)

create_igraph_from_matrix(matrix_fa_combined_1_3, save_name = "fa_combined_1_3")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K16792", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K03420", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K14098", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K07558", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K14115", taxa_level = "G")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00757")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00582")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01649")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01596")

p_plot_taxa_abundance(taxa_name = "Treponema_D")
p_plot_taxa_abundance(taxa_name = "Methanobrevibacter_A")


fa_combined_1_3_sub1 <- filter_submatrix(matrix_fa_combined_1_3,
                                         vector = c("fa_p", "hg_cp", "fa_ca",
                                                    "p_K03231", "p_K01649", "p_K00582", "p_K00123", "p_K00757",
                                                    "p_K01596", "p_K00609", "p_K14126", "p_K21990", "p_K01191",
                                                    "p_Clostridium_AI", 
                                                    "g_K03420", "g_K14098", "g_K07558", "g_K16792", "g_K14115",
                                                    "ssc_397108", "ssc_397376", "ssc_100521982", 
                                                    "A0A9D5B9N7_PEA", 
                                                    "m_acetate", "m_valerate", "m_uracil", "m_propionate", 
                                                    "m_butyrate", "m_aspartate", "m_isobutyrate"))

create_igraph_from_matrix(fa_combined_1_3_sub1, save_name = "fa_combined_1_3_sub1")

fa_combined_1_3_sub1_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_3_sub1), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_3_sub1_ko_enriched <- enrich_keggs(filter(fa_combined_1_3_sub1_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

g_plot_kegg_taxa_origin(kegg_ko_filter = "K16792", taxa_level = "S")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K03420", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K14115", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K14098", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K07558", taxa_level = "G")

g_plot_taxa_abundance(taxa_name = "Methanobacteriaceae", taxa_level = "F")
g_plot_taxa_abundance(taxa_name = "Methanobrevibacter_A", taxa_level = "G")
g_plot_taxa_abundance(taxa_name = "Methanobrevibacter_A sp900769095", taxa_level = "S")
g_plot_taxa_abundance(taxa_name = "Methanobrevibacter_B", taxa_level = "G")
g_plot_taxa_abundance(taxa_name = "Methanosphaera", taxa_level = "G")

p_plot_taxa_abundance(taxa_name = "Methanobrevibacter_A")
p_plot_taxa_abundance(taxa_name = "Faecousia")
p_plot_taxa_abundance(taxa_name = "Clostridium_AI")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K00123", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00582")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K14126", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01191")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00609")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K21990", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K00757")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01649")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01596", taxa_level = "G")

ssc_plot_kegg_abundance(kegg_ko_filter = "397108")
ssc_plot_kegg_abundance(kegg_ko_filter = "397376")
ssc_plot_kegg_abundance(kegg_ko_filter = "100521982")

p_plot_kegg_for_taxa(taxa_name = "Methanobrevibacter_A", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Clostridium_AI", use_names = T)

archaea <- breport_rel_abd_filtered %>% 
  filter_faeces() %>%
  filter(rank == "S") %>%
  filter(str_detect(name, "Methano")) %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd))

ggplot(archaea, aes(x = diet, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors)

fa_combined_1_3_sub2 <- filter_submatrix(matrix_fa_combined_1_3,
                                         vector = c("fa_cp",
                                                    "p_K18676", "p_K01915", "p_K02035", "p_K02051",
                                                    "p_PeH17",
                                                    "g_K01031", "g_K04835", "g_K01846", "g_K00772", "g_K01758",
                                                    "ssc_396921", "ssc_445518", "ssc_396674", "ssc_445532",
                                                    "A0A9D4YG61_PEA", "A0A9D4W4Y0_PEA", "A0A9D4W9W0_PEA",
                                                    "A0A9D5AY78_PEA", "A0A9D5A4B1_PEA",
                                                    "m_methionine", "m_phenylalanine", "m_tyrosine", "m_valine"))

fa_combined_1_3_sub2_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_3_sub2), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_3_sub2_ko_enriched <- enrich_keggs(filter(fa_combined_1_3_sub2_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

pea_plot_protein_abundance(protein_name = "A0A9D4YG61_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4W9W0_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4W4Y0_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5AY78_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5A4B1_PEA")

g_plot_taxa_abundance(taxa_name = "Onthomorpha", taxa_level = "G")
gtdb <- readRDS("clean/3_gtdb_v220_taxonomy.RDS")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K01758")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K04835")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K01846")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K01031", taxa_level = "G")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K00772", taxa_level = "F")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K16363")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K01267")
g_plot_kegg_taxa_origin(kegg_ko_filter = "K02535")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K01915")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02035")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02051")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K18676")

p_plot_taxa_abundance(taxa_name = "PeH17")
p_plot_taxa_abundance(taxa_name = "Phil1")

ssc_plot_kegg_abundance(kegg_ko_filter = "445532")
ssc_plot_kegg_abundance(kegg_ko_filter = "445518")
ssc_plot_kegg_abundance(kegg_ko_filter = "396674")

peh17 <- p_plot_kegg_for_taxa(taxa_name = "PeH17", use_names = T)
p_plot_kegg_for_taxa(taxa_name = "Phil1", use_names = T)

fa_combined_1_3_sub2_1 <- fa_combined_1_3_sub2[fa_combined_1_3_sub2["p_PeH17",] != 0,
                                               fa_combined_1_3_sub2[,"p_PeH17"] != 0]
create_igraph_from_matrix(fa_combined_1_3_sub2_1)

# 1 vs 4

fa_combined_1_4 <- filter_input_smccn(fa_combined, meta_fa, include = c("1", "4"))

cor_mean <- diablo_pairwise_correlations(fa_combined_1_4[[1]])

diablo_fa_combined_1_4 <- diablo_auto_model(omics_list = fa_combined_1_4[[1]], meta = fa_combined_1_4[[2]]$diet,
                                            design_value = cor_mean, ncomp_value = 2) # 2 for plotting

matrix_fa_combined_1_4 <- plot_diablo(diablo_fa_combined_1_4, cutoff = 0.8, save_name = "fa_combined_1_4")
#matrix_fa_combined_1_4 <- plot_diablo(diablo_fa_combined_1_4, return_full_matrix = T)

fa_combined_1_4_sub1 <- filter_submatrix(matrix_fa_combined_1_4,
                                         vector = c("fa_k", "hg_dm", "hg_tdf", "fa_ti", 
                                                    "p_K06410", "p_K04079", "p_K02652", "p_K02662", "p_K02243",
                                                    "g_K01267",
                                                    "g_UMGS124 sp019420325", "g_UMGS124 sp902464015", 
                                                    "g_HGM13006 sp029012465", "g_UMGS124 sp900555105",
                                                    "ssc_100525899", "ssc_407610", "ssc_100037943", "ssc_445461",
                                                    "ssc_397397",
                                                    "A0A9D5BQA3_PEA", "A0A9D5BN09_PEA", "A0A9D4XTC5_PEA", 
                                                    "A0A9D4VJF2_PEA", "A0A9D4WCI6_PEA",
                                                    "m_acetate", "m_butyrate", "m_propionate", "m_methionine",
                                                    "m_galactose"))

fa_combined_1_4_sub1_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_4_sub1), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_4_sub1_ko_enriched <- enrich_keggs(filter(fa_combined_1_4_sub1_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_taxa_abundance(taxa_name = "Ruminococcus")
p_plot_taxa_abundance(taxa_name = "Turicibacter")

p_plot_kegg_taxa_origin(kegg_ko_filter = "K02662")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02243")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K02652", taxa_level = "S")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K04079")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01267")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K06410", taxa_level = "S")

g_plot_kegg_taxa_origin(kegg_ko_filter = "K01267")

g_plot_taxa_abundance(taxa_name = "UMGS124 sp900555105")
g_plot_taxa_abundance(taxa_name = "UMGS124 sp902464015")
g_plot_taxa_abundance(taxa_name = "UMGS124 sp019420325")
g_plot_taxa_abundance(taxa_name = "HGM13006 sp029012465")
g_plot_taxa_abundance(taxa_name = "Coriobacteriales", taxa_level = "O")

pea_plot_protein_abundance(protein_name = "A0A9D5BN09_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4WCI6_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4VJF2_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D5BQA3_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4XTC5_PEA")

ssc_plot_kegg_abundance(kegg_ko_filter = "397397")
ssc_plot_kegg_abundance(kegg_ko_filter = "445461")
ssc_plot_kegg_abundance(kegg_ko_filter = "100037943")
ssc_plot_kegg_abundance(kegg_ko_filter = "407610")

p_plot_kegg_for_taxa(taxa_name = "Cryptobacteroides")
p_plot_kegg_for_taxa(taxa_name = "Treponema_D", use_names = T, threshold = .75)
p_plot_kegg_for_taxa(taxa_name = "Turicibacter", use_names = T)
p_plot_kegg_fot_taxa(taxa_name = "Ruminococcaceae", use_names = T, taxa_level = "F")

fa_combined_1_4_sub2 <- filter_submatrix(matrix_fa_combined_1_4,
                                         vector = c("enz_carb", "enz_try", 
                                                    "p_K21471", "p_K01447", "p_K13694", "p_K19224", "p_K01278",
                                                    "p_K04488",
                                                    "g_Roseburia sp902781225", "g_CALXSC01 sp934747265",
                                                    "g_Ruminiclostridium_E sp945921515", "g_UBA10677 sp934270565",
                                                    "g_CAG-273 sp003507395",
                                                    "ssc_397520", "ssc_100158015", "ssc_100737088", "ssc_397080",
                                                    "ssc_397012",
                                                    "A0A9D4WJM1_PEA", "A0A9D4X6Z1_PEA", "A0A9D4WX68_PEA",
                                                    "A0A9D4Y3P4_PEA", "A0A9D5BFI9_PEA",
                                                    "m_isovalerate", "m_phenylacetate", "m_isobutyrate"))

fa_combined_1_4_sub2_ko <- annotate_keggs(str_remove(rownames(fa_combined_1_4_sub2), "g_|p_") %>%
                                            str_replace("ssc_", "ssc:"))
fa_combined_1_4_sub2_ko_enriched <- enrich_keggs(filter(fa_combined_1_4_sub2_ko, str_detect(kegg_ko, "^K"))$kegg_ko)

p_plot_kegg_taxa_origin(kegg_ko_filter = "K21471")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01447")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K13694")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K04488")
p_plot_kegg_taxa_origin(kegg_ko_filter = "K01278")

g_plot_taxa_abundance(taxa_name = "CALXSC01 sp934747265")
g_plot_taxa_abundance(taxa_name = "Roseburia sp902781225")
g_plot_taxa_abundance(taxa_name = "Ruminiclostridium_E sp945921515")

pea_plot_protein_abundance(protein_name = "A0A9D4Y3P4_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4X6Z1_PEA")
pea_plot_protein_abundance(protein_name = "A0A9D4WJM1_PEA")

ssc_plot_kegg_abundance(kegg_ko_filter = "397520")
ssc_plot_kegg_abundance(kegg_ko_filter = "397080")
ssc_plot_kegg_abundance(kegg_ko_filter = "397012")
ssc_plot_kegg_abundance(kegg_ko_filter = "397520")
ssc_plot_kegg_abundance(kegg_ko_filter = "100158015")

# save objects
save.image("temp/74_multiomics_diablo_faeces.RData")


