library(here)
library(patchwork)
library(cowplot)

source(here("30_omics_functions.R"))

source("E:/R/source/ggplot2_theme_JAS.R")

# load meta

meta <- readRDS("clean/meta1.RDS") %>%
  mutate(sampleno = str_extract(sampleid, "[0-9]{2}$")) %>%
  mutate(diet = case_when(diet == "1" ~ "SP",
                          diet == "2" ~ "WP1",
                          diet == "3" ~ "WP2",
                          diet == "4" ~ "SFP"),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP")))

# Metagenomics plots

# Alpha diversity

# load data

bracken_reads <- readRDS("clean/3_bracken_reads_long_reads.RDS")%>% filter(sampleid != "103")
bracken_rel_abd <- readRDS("clean/3_bracken_reads_long_rel_abd.RDS")%>% filter(sampleid != "103")

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
  inner_join(table_simpson_fa, by = "diet") %>%
  pivot_longer(-diet, names_to = "index", values_to = "value") %>%
  pivot_wider(names_from = "diet", values_from = "value") %>%
  pivot_longer(c("SP", "WP1", "WP2", "SFP"), names_to = "diet", values_to = "mean") %>%
  separate(index, into = c("index", "region"), sep = "_") %>%
  mutate(matrix = ifelse(region == "il", "ileal digesta", "faeces")) %>%
  mutate(diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP")))

alpha_combined <- alpha %>%
  pivot_longer(c(shannon, observed, simpson), names_to = "index", values_to = "value") %>%
  left_join(table_alpha, by = c("matrix", "index", "diet")) %>%
  filter(index != "simpson") %>%
  mutate(index = case_when(index == "observed" ~ "Observed species",
                           index == "shannon" ~ "Shannon index")) %>%
  mutate(mean = str_extract(mean, "[a-z]+")) %>% # only letters
  group_by(matrix, index) %>%
  mutate(max_value = max(value)) %>% # for p value y value
  ungroup() %>%
  group_by(matrix, index, diet) %>%
  mutate(max_value_letter = max(value)) %>% # for letter y value
  ungroup()

a <- ggplot(filter(alpha_combined, matrix == "ileal digesta"), aes(x = diet, y = value, fill = diet)) +
  geom_boxplot(outliers = F, position = position_dodge(width = .9), show.legend = F, alpha = .5) +
  geom_quasirandom(aes(), dodge.width = .9, show.legend = F, width = .2) +
  facet_wrap(~index, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Diet") +
  geom_text(aes(x = 3.5, y = max_value*1.05, label = paste0("italic(P) == ", `P-value`)), 
            parse = TRUE, size = 6, family = "arial") +
  theme(axis.text.x=element_text(size=16, color="black"),
        strip.text = element_text(size=18, color = "black"))

b <- ggplot(filter(alpha_combined, matrix == "faeces"), aes(x = diet, y = value, fill = diet)) +
  geom_boxplot(outliers = F, position = position_dodge(width = .9), show.legend = F, alpha = .5) +
  geom_quasirandom(aes(), dodge.width = .9, show.legend = F, width = .2) +
  facet_wrap(~index, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(y = "", x = "Diet") +
  geom_text(aes(x = diet, y = max_value_letter*1.02, label = mean), 
            parse = TRUE, size = 5, family = "arial") +
  geom_text(aes(x = 3.5, y = max_value*1.05, label = paste0("italic(P) == ", `P-value`)), 
            parse = TRUE, size = 6, family = "arial") +
  theme(axis.text.x=element_text(size=16, color="black"),
        strip.text = element_text(size=18, color = "black"))

# beta diversity

c <- do_ordination(bracken_rel_abd, region = "ileum", title = "Taxonomy reads ileum", save = F) +
  labs(title = "", color = "Diet", shape = "Diet")

d <- do_ordination(bracken_rel_abd, region = "faeces", title = "Taxonomy reads ileum", save = F) +
  labs(title = "", color = "Diet", shape = "Diet") +
  scale_color_manual(values = colors, labels = c("SP", "WP1", "WP2", "SFP")) +
  scale_shape_manual(values = c(15,16,17,18,21,22,23,24), labels = c("SP", "WP1", "WP2", "SFP"))

# taxonomic composition


breport_rel_abd_filtered <- readRDS("clean/3_breport_reads_long_rel_abd_filtered.RDS")%>% filter(sampleid != "103")

taxa_barplot_from_ktable(breport_rel_abd_filtered, meta = meta, selected_rank = "G", selected_matrix = "ileal digesta", 
                         title = "", save = F)

filtered_table_il_g <- filter_ktable(ktable = breport_rel_abd_filtered, meta = meta, selected_rank = "G", 
                                selected_matrix = "ileal digesta")
aggregated_table_il_g <- aggregate_low_abundant_taxa(filtered_table = filtered_table_il_g, threshold = 1)

filtered_table_fa_g <- filter_ktable(ktable = breport_rel_abd_filtered, meta = meta, selected_rank = "G", 
                                   selected_matrix = "faeces")
aggregated_table_fa_g <- aggregate_low_abundant_taxa(filtered_table = filtered_table_fa_g, threshold = 1)

######### Metaproteomics for colors

taxonomy_norm_imp_rel_abd <- readRDS("clean/4_taxonomy_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

filtered_table_il_p <- filter_ktable(ktable = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", 
                                   selected_matrix = "ileal digesta")
aggregated_table_il_p <- aggregate_low_abundant_taxa(filtered_table = filtered_table_il_p, threshold = 1)

filtered_table_fa_p <- filter_ktable(ktable = taxonomy_norm_imp_rel_abd, meta = meta, selected_rank = "G", 
                                   selected_matrix = "faeces")
aggregated_table_fa_p <- aggregate_low_abundant_taxa(filtered_table = filtered_table_fa_p, threshold = 1)

################

color_vector <- c(sort(unique(c(aggregated_table_il_g$name, aggregated_table_fa_g$name,
                                 aggregated_table_il_p$name, aggregated_table_fa_p$name)))) 
color_df <- color_vector %>%
  as_tibble() %>%
  dplyr::rename(name = value) %>%
  filter(!name == "other") %>%
  add_row(name = "other") %>%
  add_column(color = c(colors[1:length(color_vector)-1], "grey")) %>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

plotting_table <- aggregated_table_il_g %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop")%>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

e <- ggplot(plotting_table, aes(x = diet, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "Diet", y = "Relative abundance (%)", fill = "Genera", subtitle = "") +
  scale_fill_manual(values = filter(color_df, name %in% plotting_table$name)$color) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0))

plotting_table <- aggregated_table_fa_g %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop")%>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

f <- ggplot(plotting_table, aes(x = diet, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "Diet", y = "Relative abundance (%)", fill = "Genera", subtitle = "") +
  scale_fill_manual(values = filter(color_df, name %in% plotting_table$name)$color) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0))

cd <- (c | d) + plot_layout(guides = "collect")

fig1 <- (a | b) /
  cd /
  (e | f) +
  plot_layout(heights = c(3,3,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_fig1.jpeg",
       plot = fig1,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

# differential abundance reads level

breport_reads_filtered <- readRDS("clean/3_breport_reads_long_reads_filtered.RDS") %>% filter(sampleid != "103")

# ileal digesta genus

filtered_table <- filter_ktable(ktable = breport_reads_filtered, meta = meta, selected_rank = "G", 
                                selected_matrix = "ileal digesta")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "ileal digesta") %>%
  column_to_rownames("sampleid")

# perform ancombc
output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) 

res_pair_list <- list()
comparisons <- c("dietWP1", "dietWP2", "dietSFP", "dietWP2_dietWP1", "dietSFP_dietWP1", "dietSFP_dietWP2")
comparisons_renamed <- c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")

for (i in 1:length(comparisons)) {
  current_comparison <- comparisons[i]
  current_lfc <- str_c("lfc_", current_comparison)
  current_se <- str_c("se_", current_comparison)
  current_p <- str_c("p_", current_comparison)
  current_q <- str_c("q_", current_comparison)
  current_ss <- str_c("passed_ss_", current_comparison)
  
  res_pair_filtered <- output$res_pair %>%
    dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                  q = !!sym(current_q), ss = !!sym(current_ss)) %>%
    mutate(comparison = comparisons_renamed[i]) %>%
    mutate(comparison = factor(comparison, levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  res_pair_list[[i]] <- res_pair_filtered
  names(res_pair_list)[i] <- comparisons_renamed[i]
}

p1 <- create_waterfall_plot_from_list(res_pair_list)
figS1a <- p1[[2]] + 
  theme(axis.text.y = element_text(size = 15)) | 
  p1[[3]] + 
  theme(axis.text.y = element_text(size = 16))

figS1b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  scale_size_manual(values = c(2,4,4,4,4))

legend <- tibble(Significance = c("positive LFC and passed sensitivity analysis",
                       "positive LFC and not passed sensitivity analysis",
                       "negative LFC and passed sensitivity analysis",
                       "negative LFC and not passed sensitivity analysis")) %>%
  mutate(Significant = factor(Significance, levels = c("positive LFC and passed sensitivity analysis",
                                                      "positive LFC and not passed sensitivity analysis",
                                                      "negative LFC and passed sensitivity analysis",
                                                      "negative LFC and not passed sensitivity analysis"))) %>%
  ggplot(aes(x = 1, y = 1,fill = Significance)) +
  geom_tile() +
  scale_fill_manual(values = c("#33a02c", "#b2df8a", "#e31a1c", "#fb9a99")) +
  guides(fill=guide_legend(ncol=2))
legend <- get_legend(legend)
ggdraw(legend)

figS1 <- figS1a / figS1b / legend +
  plot_layout(heights = c(3,6,1)) +
  plot_annotation(tag_levels = list(c("a", "", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS1.jpeg",
       plot = figS1,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 15,
       scale=2,
       dpi=300)

#ileal digesta species

filtered_table <- filter_ktable(ktable = breport_reads_filtered, meta = meta, selected_rank = "S", 
                                selected_matrix = "ileal digesta")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "ileal digesta") %>%
  column_to_rownames("sampleid")

# perform ancombc
output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) 

res_pair_list <- list()
comparisons <- c("dietWP1", "dietWP2", "dietSFP", "dietWP2_dietWP1", "dietSFP_dietWP1", "dietSFP_dietWP2")
comparisons_renamed <- c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")

for (i in 1:length(comparisons)) {
  current_comparison <- comparisons[i]
  current_lfc <- str_c("lfc_", current_comparison)
  current_se <- str_c("se_", current_comparison)
  current_p <- str_c("p_", current_comparison)
  current_q <- str_c("q_", current_comparison)
  current_ss <- str_c("passed_ss_", current_comparison)
  
  res_pair_filtered <- output$res_pair %>%
    dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                  q = !!sym(current_q), ss = !!sym(current_ss)) %>%
    mutate(comparison = comparisons_renamed[i]) %>%
    mutate(comparison = factor(comparison, levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  res_pair_list[[i]] <- res_pair_filtered
  names(res_pair_list)[i] <- comparisons_renamed[i]
}

p1 <- create_waterfall_plot_from_list(res_pair_list)
figS2a <- (p1[[1]] + 
  theme(axis.text.y = element_text(size = 15)) | 
  p1[[2]] + 
  theme(axis.text.y = element_text(size = 16)) |
  p1[[3]] + 
  theme(axis.text.y = element_text(size = 16))) /
  (p1[[4]] + 
  theme(axis.text.y = element_text(size = 16)) |
  p1[[5]] + 
  theme(axis.text.y = element_text(size = 16)))

figS2b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  scale_size_manual(values = c(2,4,4,4,4))


figS2 <- figS2a / figS2b / legend +
  plot_layout(heights = c(3,3,6,1)) +
  plot_annotation(tag_levels = list(c("a", "", "", "", "", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS2.jpeg",
       plot = figS2,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 20,
       scale=2,
       dpi=300)

# faeces genus level

filtered_table <- filter_ktable(ktable = breport_reads_filtered, meta = meta, selected_rank = "G", 
                                selected_matrix = "faeces")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "faeces") %>%
  column_to_rownames("sampleid")

# perform ancombc
output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) 

res_pair_list <- list()
comparisons <- c("dietWP1", "dietWP2", "dietSFP", "dietWP2_dietWP1", "dietSFP_dietWP1", "dietSFP_dietWP2")
comparisons_renamed <- c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")

for (i in 1:length(comparisons)) {
  current_comparison <- comparisons[i]
  current_lfc <- str_c("lfc_", current_comparison)
  current_se <- str_c("se_", current_comparison)
  current_p <- str_c("p_", current_comparison)
  current_q <- str_c("q_", current_comparison)
  current_ss <- str_c("passed_ss_", current_comparison)
  
  res_pair_filtered <- output$res_pair %>%
    dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                  q = !!sym(current_q), ss = !!sym(current_ss)) %>%
    mutate(comparison = comparisons_renamed[i]) %>%
    mutate(comparison = factor(comparison, levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  res_pair_list[[i]] <- res_pair_filtered
  names(res_pair_list)[i] <- comparisons_renamed[i]
}

p1 <- create_waterfall_plot_from_list(res_pair_list)
figS3a <- (p1[[1]] + 
             theme(axis.text.y = element_text(size = 16)) | 
             p1[[2]] + 
             theme(axis.text.y = element_text(size = 16)) |
             p1[[3]] + 
             theme(axis.text.y = element_text(size = 16))) /
  (p1[[4]] + 
     theme(axis.text.y = element_text(size = 16)) |
     p1[[5]] + 
     theme(axis.text.y = element_text(size = 12)) |
     p1[[6]] + 
     theme(axis.text.y = element_text(size = 16)))

figS3b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  scale_size_manual(values = c(2,4,4,4,4))


figS3 <- figS3a / figS3b / legend +
  plot_layout(heights = c(4,4,4,1)) +
  plot_annotation(tag_levels = list(c("a", "", "", "", "", "", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS3.jpeg",
       plot = figS3,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

# faeces species level

filtered_table <- filter_ktable(ktable = breport_reads_filtered, meta = meta, selected_rank = "S", 
                                selected_matrix = "faeces")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "faeces") %>%
  column_to_rownames("sampleid")

# perform ancombc
output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) 

res_pair_list <- list()
comparisons <- c("dietWP1", "dietWP2", "dietSFP", "dietWP2_dietWP1", "dietSFP_dietWP1", "dietSFP_dietWP2")
comparisons_renamed <- c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")

for (i in 1:length(comparisons)) {
  current_comparison <- comparisons[i]
  current_lfc <- str_c("lfc_", current_comparison)
  current_se <- str_c("se_", current_comparison)
  current_p <- str_c("p_", current_comparison)
  current_q <- str_c("q_", current_comparison)
  current_ss <- str_c("passed_ss_", current_comparison)
  
  res_pair_filtered <- output$res_pair %>%
    dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                  q = !!sym(current_q), ss = !!sym(current_ss)) %>%
    mutate(comparison = comparisons_renamed[i]) %>%
    mutate(comparison = factor(comparison, levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  res_pair_list[[i]] <- res_pair_filtered
  names(res_pair_list)[i] <- comparisons_renamed[i]
}

p1 <- create_waterfall_plot_from_list(res_pair_list)
figS4a <- (p1[[1]] + 
             theme(axis.text.y = element_text(size = 16)) | 
             p1[[2]] + 
             theme(axis.text.y = element_text(size = 16)) |
             p1[[3]] + 
             theme(axis.text.y = element_text(size = 16))) /
  (p1[[4]] + 
     theme(axis.text.y = element_text(size = 16)) |
     p1[[5]] + 
     theme(axis.text.y = element_text(size = 12)) |
     p1[[6]] + 
     theme(axis.text.y = element_text(size = 16)))

figS4b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  scale_size_manual(values = c(2,4,4,4,4))


figS4 <- figS4b / legend +
  plot_layout(heights = c(4,1))

ggsave(filename= "92_figS4.jpeg",
       plot = figS4,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 15,
       scale=2,
       dpi=300)

# metagenomics contig level for functions

kegg_reads_filtered <- readRDS("clean/3_kegg_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")

# filter
filtered_df <- filter_ileum(kegg_reads_filtered)

edger_object <- create_edger_object(filtered_df, abundance_column = "reads")

edger_results <- edger_analysis(edger_object = edger_object, 
                                feature_column = colnames(filtered_df)[1])

edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = colnames(filtered_df)[1])

# change contrasts in edger results
for (i in 1:length(edger_results)) {
  table <- edger_results[[i]][["table_sig"]] %>%
    mutate(contrast = case_when(contrast == "diet2-diet1" ~ "WP1 vs. SP",
                                contrast == "diet3-diet1" ~ "WP2 vs. SP",
                                contrast == "diet4-diet1" ~ "SFP vs. SP",
                                contrast == "diet3-diet2" ~ "WP2 vs. WP1",
                                contrast == "diet4-diet2" ~ "SFP vs. WP1",
                                contrast == "diet4-diet3" ~ "SFP vs. WP2"),
           contrast = factor(contrast, 
                             levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", 
                                        "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  edger_results[[i]][["table_sig"]] <- table # safe table to list for later
}

figS5a <- create_volcano_plot_from_edger(edger_results = edger_results, title = "", no_label = 5)  +
  scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
  scale_size_manual(values = c(0.5,3), labels = c("q ≥ 0.05", "q < 0.05"))


# filter faeces
filtered_df <- filter_faeces(kegg_reads_filtered)

edger_object <- create_edger_object(filtered_df, abundance_column = "reads")

edger_results <- edger_analysis(edger_object = edger_object, 
                                feature_column = colnames(filtered_df)[1])

edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = colnames(filtered_df)[1])

# change contrasts in edger results
for (i in 1:length(edger_results)) {
  table <- edger_results[[i]][["table_sig"]] %>%
    mutate(contrast = case_when(contrast == "diet2-diet1" ~ "WP1 vs. SP",
                                contrast == "diet3-diet1" ~ "WP2 vs. SP",
                                contrast == "diet4-diet1" ~ "SFP vs. SP",
                                contrast == "diet3-diet2" ~ "WP2 vs. WP1",
                                contrast == "diet4-diet2" ~ "SFP vs. WP1",
                                contrast == "diet4-diet3" ~ "SFP vs. WP2"),
           contrast = factor(contrast, 
                             levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", 
                                        "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  edger_results[[i]][["table_sig"]] <- table # safe table to list for later
}

figS5b <- create_volcano_plot_from_edger(edger_results = edger_results, title = "", no_label = 5) +
  scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
  scale_size_manual(values = c(0.5,3), labels = c("q ≥ 0.05", "q < 0.05"))



figS5 <- figS5a / figS5b +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = list(c("a", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"),
        legend.position = "bottom")

ggsave(filename= "92_figS5.jpeg",
       plot = figS5,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

####################
# Metaproteomics
###################

# micro host ratio

proteins <- readRDS("clean/4_proteins_long_raw_intensity.RDS") %>% filter(sampleid != "103")

proteins_origin <- proteins %>%
  filter(origin %in% c("pig", "micro")) %>% # only host and microbial proteins
  group_by(sampleid, origin) %>%
  summarize(intensity = sum(intensity), .groups = "drop") %>%
  group_by(sampleid) %>%
  mutate(rel_abd = intensity/sum(intensity)) %>%
  ungroup() %>%
  left_join(meta, by = "sampleid")


proteins_origin_ratio <- proteins_origin %>%
  dplyr::select(-intensity) %>%
  pivot_wider(values_from = "rel_abd", names_from = "origin") %>%
  mutate(ratio = micro/pig)


fa_ratio <- combined_comparison(select_response(filter(proteins_origin_ratio, matrix == "faeces"), "ratio"), 
                                transformation = "test") # sig lower ratio in diet4

table_ratio_fa <- create_results_table(select_response(filter(proteins_origin_ratio, matrix == "faeces"), "ratio"),
                                     fa_ratio, response = "ratio_fa")

il_ratio <- combined_comparison(select_response(filter(proteins_origin_ratio, matrix == "ileal digesta"), "ratio"), 
                                transformation = "test") # sig lower ratio in diet4

table_ratio_il <- create_results_table(select_response(filter(proteins_origin_ratio, matrix == "ileal digesta"), "ratio"),
                                       il_ratio, response = "ratio_il")

table_ratio <- table_ratio_fa %>%
  inner_join(table_ratio_il, by = "diet") %>%
  pivot_longer(-diet, names_to = "ratio", values_to = "value") %>%
  pivot_wider(names_from = "diet", values_from = "value") %>%
  pivot_longer(c("SP", "WP1", "WP2", "SFP"), names_to = "diet", values_to = "mean") %>%
  separate(ratio, into = c("ratio", "region"), sep = "_") %>%
  mutate(matrix = ifelse(region == "il", "ileal digesta", "faeces")) %>%
  dplyr::select(-ratio, -region) %>%
  mutate(diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  mutate(`P-value` = ifelse(str_detect(`P-value`, "<"), `P-value`, paste("==", `P-value`))) # == for compatibility with < and parse
 
ratio_combined <- proteins_origin_ratio %>%
  left_join(table_ratio, by = c("matrix", "diet")) %>%
  mutate(matrix = ifelse(matrix == "faeces", "Faeces", "Ileal digesta"),
         matrix = factor(matrix, levels = c("Ileal digesta", "Faeces"))) %>%
  mutate(mean = str_extract(mean, "[a-z]+")) %>% # only letters
  group_by(matrix) %>%
  mutate(max_value = max(ratio),
         min_value = min(ratio)) %>% # for p value y value
  ungroup() %>%
  group_by(matrix, diet) %>%
  mutate(max_value_letter = max(ratio)) %>% # for letter y value
  ungroup()

fig2a <- ggplot(ratio_combined, aes(x = diet, y = ratio, fill = diet)) +
  geom_boxplot(outliers = F, position = position_dodge(width = .9), show.legend = F, alpha = .5) +
  geom_quasirandom(aes(), dodge.width = .9, show.legend = F, width = .2) +
  facet_wrap(~matrix, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(y = "Microbial - host protein ratio", x = "Diet") +
  geom_text(aes(x = diet, y = max_value_letter + (max_value - min_value)*.05, label = mean), 
            parse = TRUE, size = 6, family = "arial") +
  geom_text(aes(x = 3.5, y = max_value*1.05, label = paste0("italic(P) ", `P-value`)), 
            parse = TRUE, size = 6, family = "arial") +
  theme(axis.text.x=element_text(size=16, color="black"),
        strip.text = element_text(size=18, color = "black"))

# shannon

taxa_raw_rel_abd <- readRDS("clean/4_taxa_long_raw_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

matrix <- taxa_raw_rel_abd %>%
  dplyr::select(bin, sampleid, rel_abd) %>%
  pivot_wider(names_from = "bin", values_from = "rel_abd") %>%
  column_to_rownames("sampleid") %>%
  mutate_all(~ifelse(is.na(.), 0, .))

shannon <- tibble(sampleid = rownames(matrix),
                  shannon = diversity(matrix, index = "shannon"),
                  simpson = diversity(matrix, index = "invsimpson")) %>%
  inner_join(meta, by = "sampleid")

alpha_il <- select_matrix(shannon, "ileal digesta")
alpha_fa <- select_matrix(shannon, "faeces")

comp_shannon_il <- combined_comparison(select_response(alpha_il, "shannon"), transformation = "test") # not sig
comp_shannon_fa <- combined_comparison(select_response(alpha_fa, "shannon"), transformation = "test") #not sig

table_shannon_il <- create_results_table(select_response(alpha_il, "shannon"), comp_shannon_il, response = "shannon_il")
table_shannon_fa <- create_results_table(select_response(alpha_fa, "shannon"), comp_shannon_fa, response = "shannon_fa")

table_alpha <- table_shannon_il %>%
  inner_join(table_shannon_fa, by = "diet") %>%
  pivot_longer(-diet, names_to = "index", values_to = "value") %>%
  pivot_wider(names_from = "diet", values_from = "value") %>%
  pivot_longer(c("SP", "WP1", "WP2", "SFP"), names_to = "diet", values_to = "mean") %>%
  separate(index, into = c("index", "region"), sep = "_") %>%
  mutate(matrix = ifelse(region == "il", "ileal digesta", "faeces")) %>%
  dplyr::select(-index, -region) %>%
  mutate(diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP")))

alpha_combined <- shannon %>%
  left_join(table_alpha, by = c("matrix", "diet")) %>%
  mutate(matrix = ifelse(matrix == "faeces", "Faeces", "Ileal digesta"),
         matrix = factor(matrix, levels = c("Ileal digesta", "Faeces"))) %>%
  mutate(mean = str_extract(mean, "[a-z]+")) %>% # only letters
  group_by(matrix) %>%
  mutate(max_value = max(shannon),
         min_value = min(shannon)) %>% # for p value y value
  ungroup() %>%
  group_by(matrix, diet) %>%
  mutate(max_value_letter = max(shannon)) %>% # for letter y value
  ungroup()

fig2b <- ggplot(alpha_combined, aes(x = diet, y = shannon, fill = diet)) +
  geom_boxplot(outliers = F, position = position_dodge(width = .9), show.legend = F, alpha = .5) +
  geom_quasirandom(aes(), dodge.width = .9, show.legend = F, width = .2) +
  facet_wrap(~matrix, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(y = "Shannon index", x = "Diet") +
  geom_text(aes(x = 3.5, y = max_value + (max_value-min_value)*0.1, label = paste0("italic(P) == ", `P-value`)), 
            parse = TRUE, size = 6, family = "arial") +
  theme(axis.text.x=element_text(size=16, color="black"),
        strip.text = element_text(size=18, color = "black"))

# nFR

protein_func_tax <- readRDS("clean/4_function_taxa_combined_long_raw_rel_abd_filtered.RDS")
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

nfr_il <- select_matrix(nfr, "ileal digesta")
nfr_fa <- select_matrix(nfr, "faeces")

comp_nfr_il <- combined_comparison(select_response(nfr_il, "nFR"), transformation = "test") 
comp_nfr_fa <- combined_comparison(select_response(nfr_fa, "nFR"), transformation = "test")

table_nfr_il <- create_results_table(select_response(nfr_il, "nFR"), comp_nfr_il, response = "nFR_il", digits = 2)
table_nfr_fa <- create_results_table(select_response(nfr_fa, "nFR"), comp_nfr_fa, response = "nFR_fa", digits = 2)


table_nfr <- table_nfr_il %>%
  inner_join(table_nfr_fa, by = "diet") %>%
  pivot_longer(-diet, names_to = "index", values_to = "value") %>%
  pivot_wider(names_from = "diet", values_from = "value") %>%
  pivot_longer(c("SP", "WP1", "WP2", "SFP"), names_to = "diet", values_to = "mean") %>%
  separate(index, into = c("index", "region"), sep = "_") %>%
  mutate(matrix = ifelse(region == "il", "ileal digesta", "faeces")) %>%
  dplyr::select(-index, -region) %>%
  mutate(diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP")))

nfr_combined <- nfr %>%
  left_join(table_nfr, by = c("matrix", "diet")) %>%
  mutate(matrix = ifelse(matrix == "faeces", "Faeces", "Ileal digesta"),
         matrix = factor(matrix, levels = c("Ileal digesta", "Faeces"))) %>%
  mutate(mean = str_extract(mean, "[a-z]+")) %>% # only letters
  group_by(matrix) %>%
  mutate(max_value = max(nFR),
         min_value = min(nFR)) %>% # for p value y value
  ungroup() %>%
  group_by(matrix, diet) %>%
  mutate(max_value_letter = max(nFR)) %>% # for letter y value
  ungroup()

fig2c <- ggplot(nfr_combined, aes(x = diet, y = nFR, fill = diet)) +
  geom_boxplot(outliers = F, position = position_dodge(width = .9), show.legend = F, alpha = .5) +
  geom_quasirandom(aes(), dodge.width = .9, show.legend = F, width = .2) +
  facet_wrap(~matrix, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(y = "nFR", x = "Diet") +
  geom_text(aes(x = 3.5, y = max_value + (max_value-min_value)*0.1, label = paste0("italic(P) == ", `P-value`)), 
            parse = TRUE, size = 6, family = "arial") +
  theme(axis.text.x=element_text(size=16, color="black"),
        strip.text = element_text(size=18, color = "black"))

# beta diversity on proteingroup level for micro

proteins_norm_imp <- readRDS("clean/4_proteins_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

proteins_norm_imp_micro <- proteins_norm_imp %>%
  filter(origin == "micro") %>%
  recalculate_rel_abd()

# micro proteins ileum

bray <- create_distance_matrix(input_df = proteins_norm_imp_micro, region = "ileum", dissimilarity_index = "bray",
                               abundance_column = "rel_abd")

permutation_matrix <- create_permutation_matrix(bray = bray, meta = meta, region = "ileum", 
                                                treatment_col = "diet", design = "double_latin_square")

permanova <- do_permanova(bray, meta = meta, permutation_matrix = permutation_matrix)

pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pco1 <- paste0("PCo1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100,1), "%)")
pco2 <- paste0("PCo2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,1), "%)")

x <- as_tibble(pcoa$points, rownames = "sampleid") %>%
  inner_join(meta, by = "sampleid")
# calculate centroids
x_centers <- x %>%
  group_by(matrix, diet) %>%
  summarise(c1 = mean(V1), c2 = mean(V2), .groups = "drop")

Print <- ifelse(permanova$`Pr(>F)`[1] < 0.001, paste("< 0.001"), 
                paste("=", format(round(permanova$`Pr(>F)`[1],3), nsmall = 3)))
R2 <- format(round(permanova$R2[1], 3),nsmall = 3)

fig2d <- ggplot(x) +
  geom_point(aes(x = V1, y = V2, shape = diet, color = diet), size = 4) +
  geom_point(data = x_centers, mapping = aes(x = c1, y = c2,
                                             #color = diet, 
                                             shape = diet, 
                                             fill = matrix),
             show.legend = FALSE, size = 8, alpha = 0.5) +
  stat_ellipse(aes(x = V1, y = V2, color = diet),
               level = 0.9, lwd = 1.5, alpha = 0.3, show.legend = F) +
  labs(x = pco1, y = pco2, color = "Diet", shape = "Diet",
       subtitle = bquote("PERMANOVA:"~R^2 == .(R2)*"," ~ italic(P)~.(Print))) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(15,16,17,18,21,22,23,24))

# micro proteins faeces

bray <- create_distance_matrix(input_df = proteins_norm_imp_micro, region = "faeces", dissimilarity_index = "bray",
                               abundance_column = "rel_abd")

permutation_matrix <- create_permutation_matrix(bray = bray, meta = meta, region = "faeces", 
                                                treatment_col = "diet", design = "double_latin_square")

permanova <- do_permanova(bray, meta = meta, permutation_matrix = permutation_matrix)


pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pco1 <- paste0("PCo1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100,1), "%)")
pco2 <- paste0("PCo2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,1), "%)")

x <- as_tibble(pcoa$points, rownames = "sampleid") %>%
  inner_join(meta, by = "sampleid") 
# calculate centroids
x_centers <- x %>%
  group_by(matrix, diet) %>%
  summarise(c1 = mean(V1), c2 = mean(V2), .groups = "drop")

Print <- ifelse(permanova$`Pr(>F)`[1] < 0.001, paste("< 0.001"), 
                paste("=", format(round(permanova$`Pr(>F)`[1],3), nsmall = 3)))
R2 <- format(round(permanova$R2[1], 3),nsmall = 3)

fig2e <- ggplot(x) +
  geom_point(aes(x = V1, y = V2, shape = diet, color = diet), size = 4) +
  geom_point(data = x_centers, mapping = aes(x = c1, y = c2,
                                             #color = diet, 
                                             shape = diet, 
                                             fill = matrix),
             show.legend = FALSE, size = 8, alpha = 0.5) +
  stat_ellipse(aes(x = V1, y = V2, color = diet),
               level = 0.9, lwd = 1.5, alpha = 0.3, show.legend = F) +
  labs(x = pco1, y = pco2, color = "Diet", shape = "Diet",
       subtitle = bquote("PERMANOVA:"~R^2 == .(R2)*"," ~ italic(P)~.(Print))) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(15,16,17,18,21,22,23,24))

# taxonomy barplots

  
plotting_table <- aggregated_table_il_p %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop")%>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

fig2f <- ggplot(plotting_table, aes(x = diet, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "Diet", y = "Relative abundance (%)", fill = "Genera", subtitle = "") +
  scale_fill_manual(values = filter(color_df, name %in% plotting_table$name)$color) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0))

plotting_table <- aggregated_table_fa_p %>%
  inner_join(meta, by = "sampleid") %>%
  group_by(diet, name) %>%
  summarise(rel_abd = mean(rel_abd), .groups = "drop")%>%
  mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other")))

fig2g <- ggplot(plotting_table, aes(x = diet, y = rel_abd, fill = name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  labs(x = "Diet", y = "Relative abundance (%)", fill = "Genera", subtitle = "") +
  scale_fill_manual(values = filter(color_df, name %in% plotting_table$name)$color) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0,0))

# fig2

fig2de <- (fig2d | fig2e) +
  plot_layout(guides = "collect")

fig2 <- (fig2a | fig2c) /
  fig2de /
  (fig2f | fig2g) +
  plot_layout(heights = c(3,3,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_fig2.jpeg",
       plot = fig2,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

# beta diversity on proteingroup level for host

proteins_norm_imp <- readRDS("clean/4_proteins_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

proteins_norm_imp_host <- proteins_norm_imp %>%
  filter(origin == "pig") %>%
  recalculate_rel_abd()

# host proteins ileum

bray <- create_distance_matrix(input_df = proteins_norm_imp_host, region = "ileum", dissimilarity_index = "bray",
                               abundance_column = "rel_abd")

permutation_matrix <- create_permutation_matrix(bray = bray, meta = meta, region = "ileum", 
                                                treatment_col = "diet", design = "double_latin_square")

permanova <- do_permanova(bray, meta = meta, permutation_matrix = permutation_matrix)

cld <- pairwise_permanova(bray, meta = meta, region = "ileum", design = "double_latin_square")
# replace diet groups with cld's for plotting
meta_cld <- meta %>%
  inner_join(cld, by = c("diet"="treatments")) %>%
  mutate(diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP")))

pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pco1 <- paste0("PCo1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100,1), "%)")
pco2 <- paste0("PCo2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,1), "%)")

x <- as_tibble(pcoa$points, rownames = "sampleid") %>%
  inner_join(meta_cld, by = "sampleid") %>%
  mutate(treatments_cld = str_extract(treatments_cld, "[a-z]+"))

# calculate centroids
x_centers <- x %>%
  group_by(matrix, diet, treatments_cld) %>%
  summarise(c1 = mean(V1), c2 = mean(V2), .groups = "drop")

Print <- ifelse(permanova$`Pr(>F)`[1] < 0.001, paste("< 0.001"), 
                paste("=", format(round(permanova$`Pr(>F)`[1],3), nsmall = 3)))
R2 <- format(round(permanova$R2[1], 3),nsmall = 3)

fig3a <- ggplot(x) +
  geom_point(aes(x = V1, y = V2, shape = diet, color = diet), size = 4) +
  geom_point(data = x_centers, mapping = aes(x = c1, y = c2,
                                             #color = diet, 
                                             shape = diet, 
                                             fill = matrix),
             show.legend = FALSE, size = 10, alpha = 0.5) +
  geom_text(data = x_centers, mapping = aes(x = c1, y = c2, label = treatments_cld),
            color = "white", family = "arial", size = 6) +
  stat_ellipse(aes(x = V1, y = V2, color = diet),
               level = 0.9, lwd = 1.5, alpha = 0.3, show.legend = F) +
  labs(x = pco1, y = pco2, color = "Diet", shape = "Diet",
       subtitle = bquote("PERMANOVA:"~R^2 == .(R2)*"," ~ italic(P)~.(Print))) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(15,16,17,18,21,22,23,24))

# host proteins faeces

bray <- create_distance_matrix(input_df = proteins_norm_imp_host, region = "faeces", dissimilarity_index = "bray",
                               abundance_column = "rel_abd")

permutation_matrix <- create_permutation_matrix(bray = bray, meta = meta, region = "faeces", 
                                                treatment_col = "diet", design = "double_latin_square")

permanova <- do_permanova(bray, meta = meta, permutation_matrix = permutation_matrix)

cld <- pairwise_permanova(bray, meta = meta, region = "faeces", design = "double_latin_square")
# replace diet groups with cld's for plotting
meta_cld <- meta %>%
  inner_join(cld, by = c("diet"="treatments")) %>%
  mutate(diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP")))

pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pco1 <- paste0("PCo1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100,1), "%)")
pco2 <- paste0("PCo2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,1), "%)")

x <- as_tibble(pcoa$points, rownames = "sampleid") %>%
  inner_join(meta_cld, by = "sampleid") %>%
  mutate(treatments_cld = str_extract(treatments_cld, "[a-z]+"))

# calculate centroids
x_centers <- x %>%
  group_by(matrix, diet, treatments_cld) %>%
  summarise(c1 = mean(V1), c2 = mean(V2), .groups = "drop")

Print <- ifelse(permanova$`Pr(>F)`[1] < 0.001, paste("< 0.001"), 
                paste("=", format(round(permanova$`Pr(>F)`[1],3), nsmall = 3)))
R2 <- format(round(permanova$R2[1], 3),nsmall = 3)

fig3b <- ggplot(x) +
  geom_point(aes(x = V1, y = V2, shape = diet, color = diet), size = 4) +
  geom_point(data = x_centers, mapping = aes(x = c1, y = c2,
                                             #color = diet, 
                                             shape = diet, 
                                             fill = matrix),
             show.legend = FALSE, size = 10, alpha = 0.5) +
  geom_text(data = x_centers, mapping = aes(x = c1, y = c2, label = treatments_cld),
            color = "white", family = "arial", size = 6) +
  stat_ellipse(aes(x = V1, y = V2, color = diet),
               level = 0.9, lwd = 1.5, alpha = 0.3, show.legend = F) +
  labs(x = pco1, y = pco2, color = "Diet", shape = "Diet",
       subtitle = bquote("PERMANOVA:"~R^2 == .(R2)*"," ~ italic(P)~.(Print))) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(15,16,17,18,21,22,23,24))

# differential abundance host proteins

host_kegg_norm_imp <- readRDS("clean/4_host_kegg_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
host_kegg_norm_imp_rel_abd <- readRDS("clean/4_host_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

# ileum


filtered_df <- filter_ileum(host_kegg_norm_imp)

edger_object <- create_edger_object(filtered_df, abundance_column = "intensity")

edger_results <- edger_analysis(edger_object = edger_object, 
                                feature_column = "kegg_ko")

edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = "kegg_ko")

# change contrasts in edger results
for (i in 1:length(edger_results)) {
  table <- edger_results[[i]][["table_sig"]] %>%
    mutate(contrast = case_when(contrast == "diet2-diet1" ~ "WP1 vs. SP",
                                contrast == "diet3-diet1" ~ "WP2 vs. SP",
                                contrast == "diet4-diet1" ~ "SFP vs. SP",
                                contrast == "diet3-diet2" ~ "WP2 vs. WP1",
                                contrast == "diet4-diet2" ~ "SFP vs. WP1",
                                contrast == "diet4-diet3" ~ "SFP vs. WP2"),
           contrast = factor(contrast, 
                             levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", 
                                        "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  edger_results[[i]][["table_sig"]] <- table # safe table to list for later
}

fig3c <- create_volcano_plot_from_edger(edger_results = edger_results, title = "", no_label = 5) +
  scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
  scale_size_manual(values = c(0.5,3), labels = c("q ≥ 0.05", "q < 0.05"))

# faeces

filtered_df <- filter_faeces(host_kegg_norm_imp)

edger_object <- create_edger_object(filtered_df, abundance_column = "intensity")

edger_results <- edger_analysis(edger_object = edger_object, 
                                feature_column = "kegg_ko")

edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = "kegg_ko")

# change contrasts in edger results
for (i in 1:length(edger_results)) {
  table <- edger_results[[i]][["table_sig"]] %>%
    mutate(contrast = case_when(contrast == "diet2-diet1" ~ "WP1 vs. SP",
                                contrast == "diet3-diet1" ~ "WP2 vs. SP",
                                contrast == "diet4-diet1" ~ "SFP vs. SP",
                                contrast == "diet3-diet2" ~ "WP2 vs. WP1",
                                contrast == "diet4-diet2" ~ "SFP vs. WP1",
                                contrast == "diet4-diet3" ~ "SFP vs. WP2"),
           contrast = factor(contrast, 
                             levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", 
                                        "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  edger_results[[i]][["table_sig"]] <- table # safe table to list for later
}

fig3d <- create_volcano_plot_from_edger(edger_results = edger_results, title = "", no_label = 5) +
  scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
  scale_size_manual(values = c(0.5,3), labels = c("q ≥ 0.05", "q < 0.05"))

# fig3

fig3ab <- (fig3a | fig3b) +
  plot_layout(guides = "collect")

fig3cd <- (fig3c / fig3d) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig3 <- fig3ab /
  fig3cd +
  plot_layout(heights = c(4,6)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_fig3.jpeg",
       plot = fig3,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

# Ancombc metaproteomics

# ileal digesta genus

taxonomy_norm_imp_intensity <- readRDS("clean/4_taxonomy_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")

filtered_table <- filter_ktable(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "G", 
                                selected_matrix = "ileal digesta")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "ileal digesta") %>%
  column_to_rownames("sampleid")

# perform ancombc
output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) 

res_pair_list <- list()
comparisons <- c("dietWP1", "dietWP2", "dietSFP", "dietWP2_dietWP1", "dietSFP_dietWP1", "dietSFP_dietWP2")
comparisons_renamed <- c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")

for (i in 1:length(comparisons)) {
  current_comparison <- comparisons[i]
  current_lfc <- str_c("lfc_", current_comparison)
  current_se <- str_c("se_", current_comparison)
  current_p <- str_c("p_", current_comparison)
  current_q <- str_c("q_", current_comparison)
  current_ss <- str_c("passed_ss_", current_comparison)
  
  res_pair_filtered <- output$res_pair %>%
    dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                  q = !!sym(current_q), ss = !!sym(current_ss)) %>%
    mutate(comparison = comparisons_renamed[i]) %>%
    mutate(comparison = factor(comparison, levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  res_pair_list[[i]] <- res_pair_filtered
  names(res_pair_list)[i] <- comparisons_renamed[i]
}

p1 <- create_waterfall_plot_from_list(res_pair_list)
figS6a <- p1[[3]] + 
  theme(axis.text.y = element_text(size = 16)) | 
  p1[[5]] + 
  theme(axis.text.y = element_text(size = 16))

figS6b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  scale_size_manual(values = c(2,4,4,4,4))

legend <- tibble(Significance = c("positive LFC and passed sensitivity analysis",
                                  "positive LFC and not passed sensitivity analysis",
                                  "negative LFC and passed sensitivity analysis",
                                  "negative LFC and not passed sensitivity analysis")) %>%
  mutate(Significant = factor(Significance, levels = c("positive LFC and passed sensitivity analysis",
                                                       "positive LFC and not passed sensitivity analysis",
                                                       "negative LFC and passed sensitivity analysis",
                                                       "negative LFC and not passed sensitivity analysis"))) %>%
  ggplot(aes(x = 1, y = 1,fill = Significance)) +
  geom_tile() +
  scale_fill_manual(values = c("#33a02c", "#b2df8a", "#e31a1c", "#fb9a99")) +
  guides(fill=guide_legend(ncol=2))
legend <- get_legend(legend)
ggdraw(legend)

figS6 <- figS6a / figS6b / legend +
  plot_layout(heights = c(3,6,1)) +
  plot_annotation(tag_levels = list(c("a", "", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS6.jpeg",
       plot = figS6,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 15,
       scale=2,
       dpi=300)

# faeces genus

filtered_table <- filter_ktable(ktable = taxonomy_norm_imp_intensity, meta = meta, selected_rank = "G", 
                                selected_matrix = "faeces")

matrix <- prepare_table_for_ancombc(filtered_table)

filtered_meta <- filter_meta(meta = meta, selected_matrix = "faeces") %>%
  column_to_rownames("sampleid")

# perform ancombc
output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                   p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                   taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) 

res_pair_list <- list()
comparisons <- c("dietWP1", "dietWP2", "dietSFP", "dietWP2_dietWP1", "dietSFP_dietWP1", "dietSFP_dietWP2")
comparisons_renamed <- c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")

for (i in 1:length(comparisons)) {
  current_comparison <- comparisons[i]
  current_lfc <- str_c("lfc_", current_comparison)
  current_se <- str_c("se_", current_comparison)
  current_p <- str_c("p_", current_comparison)
  current_q <- str_c("q_", current_comparison)
  current_ss <- str_c("passed_ss_", current_comparison)
  
  res_pair_filtered <- output$res_pair %>%
    dplyr::select(taxon, lfc = !!sym(current_lfc), se = !!sym(current_se), p = !!sym(current_p),
                  q = !!sym(current_q), ss = !!sym(current_ss)) %>%
    mutate(comparison = comparisons_renamed[i]) %>%
    mutate(comparison = factor(comparison, levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  res_pair_list[[i]] <- res_pair_filtered
  names(res_pair_list)[i] <- comparisons_renamed[i]
}

p1 <- create_waterfall_plot_from_list(res_pair_list)
figS7a <- (p1[[1]] + 
  theme(axis.text.y = element_text(size = 16)) | 
  p1[[2]] + 
  theme(axis.text.y = element_text(size = 16)) | 
  p1[[3]] + 
  theme(axis.text.y = element_text(size = 12))) /
  (p1[[4]] + 
  theme(axis.text.y = element_text(size = 16)) | 
  p1[[5]] + 
  theme(axis.text.y = element_text(size = 14)) | 
  p1[[6]] + 
  theme(axis.text.y = element_text(size = 16)))

figS7b <- create_volcano_plot_from_list(res_pair_list, title = "") +
  scale_size_manual(values = c(2,4,4,4,4))

legend <- tibble(Significance = c("positive LFC and passed sensitivity analysis",
                                  "positive LFC and not passed sensitivity analysis",
                                  "negative LFC and passed sensitivity analysis",
                                  "negative LFC and not passed sensitivity analysis")) %>%
  mutate(Significant = factor(Significance, levels = c("positive LFC and passed sensitivity analysis",
                                                       "positive LFC and not passed sensitivity analysis",
                                                       "negative LFC and passed sensitivity analysis",
                                                       "negative LFC and not passed sensitivity analysis"))) %>%
  ggplot(aes(x = 1, y = 1,fill = Significance)) +
  geom_tile() +
  scale_fill_manual(values = c("#33a02c", "#b2df8a", "#e31a1c", "#fb9a99")) +
  guides(fill=guide_legend(ncol=2))
legend <- get_legend(legend)
ggdraw(legend)

figS7 <- figS7a / figS7b / legend +
  plot_layout(heights = c(4,4,6,1)) +
  plot_annotation(tag_levels = list(c("a", "", "","","","", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS7.jpeg",
       plot = figS7,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

# differential abundance microbial proteins

kegg_norm_imp <- readRDS("clean/4_kegg_long_norm_imp_intensity_filtered.RDS") %>% filter(sampleid != "103")
kegg_norm_imp_rel_abd <- readRDS("clean/4_kegg_long_norm_imp_rel_abd_filtered.RDS") %>% filter(sampleid != "103")

# ileum


filtered_df <- filter_ileum(kegg_norm_imp)

edger_object <- create_edger_object(filtered_df, abundance_column = "intensity")

edger_results <- edger_analysis(edger_object = edger_object, 
                                feature_column = "kegg_ko")

edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = "kegg_ko")

# change contrasts in edger results
for (i in 1:length(edger_results)) {
  table <- edger_results[[i]][["table_sig"]] %>%
    mutate(contrast = case_when(contrast == "diet2-diet1" ~ "WP1 vs. SP",
                                contrast == "diet3-diet1" ~ "WP2 vs. SP",
                                contrast == "diet4-diet1" ~ "SFP vs. SP",
                                contrast == "diet3-diet2" ~ "WP2 vs. WP1",
                                contrast == "diet4-diet2" ~ "SFP vs. WP1",
                                contrast == "diet4-diet3" ~ "SFP vs. WP2"),
           contrast = factor(contrast, 
                             levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", 
                                        "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  edger_results[[i]][["table_sig"]] <- table # safe table to list for later
}

figS8a <- create_volcano_plot_from_edger(edger_results = edger_results, title = "", no_label = 5) +
  scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
  scale_size_manual(values = c(0.5,3), labels = c("q ≥ 0.05", "q < 0.05"))

# faeces

filtered_df <- filter_faeces(kegg_norm_imp)

edger_object <- create_edger_object(filtered_df, abundance_column = "intensity")

edger_results <- edger_analysis(edger_object = edger_object, 
                                feature_column = "kegg_ko")

edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = "kegg_ko")

# change contrasts in edger results
for (i in 1:length(edger_results)) {
  table <- edger_results[[i]][["table_sig"]] %>%
    mutate(contrast = case_when(contrast == "diet2-diet1" ~ "WP1 vs. SP",
                                contrast == "diet3-diet1" ~ "WP2 vs. SP",
                                contrast == "diet4-diet1" ~ "SFP vs. SP",
                                contrast == "diet3-diet2" ~ "WP2 vs. WP1",
                                contrast == "diet4-diet2" ~ "SFP vs. WP1",
                                contrast == "diet4-diet3" ~ "SFP vs. WP2"),
           contrast = factor(contrast, 
                             levels = c("WP1 vs. SP", "WP2 vs. SP", "SFP vs. SP", 
                                        "WP2 vs. WP1", "SFP vs. WP1", "SFP vs. WP2")))
  edger_results[[i]][["table_sig"]] <- table # safe table to list for later
}

figS8b <- create_volcano_plot_from_edger(edger_results = edger_results, title = "", no_label = 5) +
  scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
  scale_size_manual(values = c(0.5,3), labels = c("q ≥ 0.05", "q < 0.05"))

figS8 <- figS8a / figS8b +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = list(c("a", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"),
        legend.position = "bottom")

ggsave(filename= "92_figS8.jpeg",
       plot = figS8,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

######
## Multiomics
###

# load image
load("temp/73_multiomics_diablo_ileum.RData")
load("temp/74_multiomics_diablo_faeces.RData")
source(here("70_multiomics_functions.R"))

create_variate_plot <- function(diablo_object, shapes) {
  variates <- diablo_object$variates
  
  for (i in 1:(length(variates)-1)) {
    variates_extract <- variates[[i]] %>%
      as_tibble(rownames = "sampleno") %>%
      mutate(block = names(variates)[i])
    if (i == 1) {
      variates_df <- variates_extract
    } else {
      variates_df <- rbind(variates_df, variates_extract)
    }
  }
  
  variates_df <- variates_df %>%
    inner_join(dplyr::select(meta, sampleno, diet) %>% distinct(), by = "sampleno") %>%
    mutate(block = factor(block, 
                          levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")))
  
  p <- ggplot(variates_df, aes(x = comp1, y = comp2, color = block, shape = diet)) +
    geom_point(size = 5, show.legend = F) +
    scale_shape_manual(values = shapes) +
    scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
    #scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")) +
    labs(x = "Component 1", y = "Component 2", color = "Block", shape = "Diet")
  return(p)
}

create_loading_plots <- function(diablo_object, shapes1, shapes2) {
  loadings <- diablo_object$loadings
  
  for (i in 1:(length(loadings)-1)) {
    loadings_extract <- loadings[[i]] %>%
      as_tibble(rownames = "sampleno") %>%
      mutate(block = names(loadings)[i])
    if (i == 1) {
      loadings_df <- loadings_extract
    } else {
      loadings_df <- rbind(loadings_df, loadings_extract)
    }
  }
  
  loadings_comp1 <- loadings_df %>%
    filter(comp1 != 0) %>%
    group_by(block) %>%
    mutate(rank = row_number(desc(abs(comp1)))) %>% # filter top 5 loadings for each block
    ungroup() %>%
    filter(rank <= 5) %>%
    mutate(rank = row_number(comp1)) %>%
    mutate(label_neg = ifelse(comp1 < 0, sampleno, ""),
           label_pos = ifelse(comp1 > 0, sampleno, "")) %>%
    mutate(block = factor(block, 
                          levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"))) %>%
    mutate(shape = ifelse(comp1 > 0, "1", "2")) # dummy shape variables
  
  a <- ggplot(loadings_comp1, aes(x = comp1, y = rank)) +
    geom_segment(aes(xend = 0), size = 1) +
    geom_point(aes(color = block, shape = shape), size = 5, show.legend = F) +
    geom_text(aes(label = label_neg, x = 0.02), hjust = "left", size = 4, family = "arial") +
    geom_text(aes(label = label_pos, x = -0.02), hjust = "right", size = 4, family = "arial") +
    scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
    scale_shape_manual(values = shapes1) +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Loadings Component 1")
  
  loadings_comp2 <- loadings_df %>%
    filter(comp2 != 0) %>%
    group_by(block) %>%
    mutate(rank = row_number(desc(abs(comp2)))) %>% # filter top 5 loadings for each block
    ungroup() %>%
    filter(rank <= 5) %>%
    mutate(rank = row_number(comp2)) %>%
    mutate(label_neg = ifelse(comp2 < 0, sampleno, ""),
           label_pos = ifelse(comp2 > 0, sampleno, "")) %>%
    mutate(block = factor(block, 
                          levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"))) %>%
    mutate(shape = ifelse(comp2 > 0, "1", "2")) # dummy shape variables
  
  b <- ggplot(loadings_comp2, aes(x = comp2, y = rank)) +
    geom_segment(aes(xend = 0), size = 1) +
    geom_point(aes(color = block, shape = shape), size = 5, show.legend = F) +
    geom_text(aes(label = label_neg, x = 0.02), hjust = "left", size = 4, family = "arial") +
    geom_text(aes(label = label_pos, x = -0.02), hjust = "right", size = 4, family = "arial") +
    scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
    scale_shape_manual(values = shapes2) +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Loadings Component 2")
  
  p <- a | b
  
  return(p)
}

# fig 9

figS9a <- create_variate_plot(diablo_il_combined_1_2, shapes = c(15,16))
figS9bc <- create_loading_plots(diablo_il_combined_1_2, shapes1 = c(15,16), shapes2 = c(15,16))

legend <- tibble(block = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"),
                 diet = rep(c("SP", "WP1"), times = 3)) %>%
  mutate(block = factor(block, 
                        levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  ggplot(aes(x = 1, y = 1,color = block, shape = diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
  #guides(fill=guide_legend(ncol=2)) +
  labs(shape = "Diet", color = "Block") +
  theme(legend.position = "bottom")
legend <- get_legend(legend)
ggdraw(legend)

figS9 <- (((figS9a | figS9bc) + plot_layout(widths = c(4,6))) / 
            (legend)) +
  plot_layout(heights = c(4,1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c","", "d", "e", "f", ""))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS9.jpeg",
       plot = figS9,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 13,
       scale=2,
       dpi=300)

# fig 10

figS10a <- create_variate_plot(diablo_il_combined_1_3, shapes = c(15,17))
figS10bc <- create_loading_plots(diablo_il_combined_1_3, shapes1 = c(15,17), shapes2 = c(15,17))

legend <- tibble(block = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"),
                 diet = rep(c("SP", "WP2"), times = 3)) %>%
  mutate(block = factor(block, 
                        levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  ggplot(aes(x = 1, y = 1,color = block, shape = diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(15,17)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
  #guides(fill=guide_legend(ncol=2)) +
  labs(shape = "Diet", color = "Block") +
  theme(legend.position = "bottom")
legend <- get_legend(legend)
ggdraw(legend)

empty <- ggplot() +
  theme_void()

figS10 <- (((figS10a | figS10bc) + plot_layout(widths = c(4,6))) / 
            (legend) / 
            empty) +
  plot_layout(heights = c(4,1,.1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c","", "d", "e", "f", ""))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS10.jpeg",
       plot = figS10,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 13,
       scale=2,
       dpi=300)


# fig 11

figS11a <- create_variate_plot(diablo_il_combined_1_4, shapes = c(15,18))
figS11bc <- create_loading_plots(diablo_il_combined_1_4, shapes1 = c(15,18), shapes2 = c(18,15))

legend <- tibble(block = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"),
                 diet = rep(c("SP", "SFP"), times = 3)) %>%
  mutate(block = factor(block, 
                        levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  ggplot(aes(x = 1, y = 1,color = block, shape = diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
  #guides(fill=guide_legend(ncol=2)) +
  labs(shape = "Diet", color = "Block") +
  theme(legend.position = "bottom")
legend <- get_legend(legend)
ggdraw(legend)

empty <- ggplot() +
  theme_void()

figS11 <- (((figS11a | figS11bc) + plot_layout(widths = c(4,6))) / 
             (legend) / 
             empty) +
  plot_layout(heights = c(4,1,.1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c","", "d", "e", "f", ""))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS11.jpeg",
       plot = figS11,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 13,
       scale=2,
       dpi=300)

# fig S12

figS12a <- create_variate_plot(diablo_fa_combined_1_2, shapes = c(15,16))
figS12bc <- create_loading_plots(diablo_fa_combined_1_2, shapes1 = c(15,16), shapes2 = c(15,16))

legend <- tibble(block = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"),
                 diet = rep(c("SP", "WP1"), times = 3)) %>%
  mutate(block = factor(block, 
                        levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  ggplot(aes(x = 1, y = 1,color = block, shape = diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(15,16)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
  #guides(fill=guide_legend(ncol=2)) +
  labs(shape = "Diet", color = "Block") +
  theme(legend.position = "bottom")
legend <- get_legend(legend)
ggdraw(legend)

empty <- ggplot() +
  theme_void()

figS12 <- (((figS12a | figS12bc) + plot_layout(widths = c(4,6))) / 
             (legend) / 
             empty) +
  plot_layout(heights = c(4,1,.1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c","", "d", "e", "f", ""))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS12.jpeg",
       plot = figS12,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 13,
       scale=2,
       dpi=300)

# fig S13

figS13a <- create_variate_plot(diablo_fa_combined_1_3, shapes = c(15,17))
figS13bc <- create_loading_plots(diablo_fa_combined_1_3, shapes1 = c(17,15), shapes2 = c(17,15))

legend <- tibble(block = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"),
                 diet = rep(c("SP", "WP2"), times = 3)) %>%
  mutate(block = factor(block, 
                        levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  ggplot(aes(x = 1, y = 1,color = block, shape = diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(15,17)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
  #guides(fill=guide_legend(ncol=2)) +
  labs(shape = "Diet", color = "Block") +
  theme(legend.position = "bottom")
legend <- get_legend(legend)
ggdraw(legend)

empty <- ggplot() +
  theme_void()

figS13 <- (((figS13a | figS13bc) + plot_layout(widths = c(4,6))) / 
             (legend) / 
             empty) +
  plot_layout(heights = c(4,1,.1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c","", "d", "e", "f", ""))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS13.jpeg",
       plot = figS13,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 13,
       scale=2,
       dpi=300)

# fig S14

figS14a <- create_variate_plot(diablo_fa_combined_1_4, shapes = c(15,18))
figS14bc <- create_loading_plots(diablo_fa_combined_1_4, shapes1 = c(15,18), shapes2 = c(18,15))

legend <- tibble(block = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics"),
                 diet = rep(c("SP", "SFP"), times = 3)) %>%
  mutate(block = factor(block, 
                        levels = c("Nutrition", "Metagenomics", "Metaproteomics", "Host", "Pea", "Metabolomics")),
         diet = factor(diet, levels = c("SP", "WP1", "WP2", "SFP"))) %>%
  ggplot(aes(x = 1, y = 1,color = block, shape = diet)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(15,18)) +
  scale_color_manual(values = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#ffff99")) +
  #guides(fill=guide_legend(ncol=2)) +
  labs(shape = "Diet", color = "Block") +
  theme(legend.position = "bottom")
legend <- get_legend(legend)
ggdraw(legend)

figS14 <- (((figS14a | figS14bc) + plot_layout(widths = c(4,6))) / 
             (legend)) +
  plot_layout(heights = c(4,1)) +
  plot_annotation(tag_levels = list(c("a", "b", "c","", "d", "e", "f", ""))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS14.jpeg",
       plot = figS14,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 13,
       scale=2,
       dpi=300)

# networks
library(qgraph)


# S10d

matrix_il_combined_1_3 <- plot_diablo(diablo_il_combined_1_3, cutoff_value = 0.85)

g <- graph_from_adjacency_matrix(matrix_il_combined_1_3, mode = "undirected", weighted = T, diag = F)
cluster <- list("1" = c("enz_chy", "pcd_gly", "pcd_arg", "pcd_ala", "pcd_lys", "pcd_asp",
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
                        "m_phenylalanine"),
                "2" = c("il_total_starch",
                        "p_K06305", "p_K03640", "p_K01176", "p_K02217",
                        "p_Basfia_A",
                        "g_K10542", "g_K03775", "g_K14441", "g_K03592", "g_K04085",
                        "ssc_397060", "ssc_100511536", "ssc_808504", "ssc_397113",
                        "A0A9D4VWB7_PEA", "A0A9D4VWB7_PEA", "A0A9D4ZWX1_PEA",
                        "A0A9D5BJS9_PEA", "A0A9D4WDI6_PEA",
                        "m_trimethylamine"))
layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
                                               vcount = vcount(g), area = vcount(g)^3, repulse.rad=vcount(g)^3.5)
layout_fr <- norm_coords(layout_fr)
# assign colors to weight values
color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
E(g)$color <- color.jet(100)[color_index]
E(g)$width <- .75
node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                          str_detect(V(g)$name, "^g_") ~ 2,
                          str_detect(V(g)$name, "^p_") ~ 3,
                          str_detect(V(g)$name, "^ssc_") ~ 4,
                          str_detect(V(g)$name, "_PEA$") ~ 5,
                          str_detect(V(g)$name, "^m_") ~ 6,
                          .default = 0)
V(g)$color <- node_color[color_vector]
V(g)$label.family <- "arial"
V(g)$label.font <- 1
V(g)$label.color <- "black"
V(g)$label.cex <- .75
jpeg(filename = paste0("plots/", get_script_number(), "_figS10d", save_name, ".jpeg"),
     width = 20, height = 15, unit="cm", res = 500, pointsize = 12, family = "arial")
par(mar = c(0,0,0,0))
plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = F,
     vertex.shape = "rectangle", 
     vertex.size = str_width(V(g)$name)*1.8,
     vertex.size2 = 3,
     vertex.frame.width = .1,
     mark.groups = cluster,
     mark.shape = .5,
     mark.expand = 0.5, 
     mark.col = NA,
     mark.border = "grey10",
     mark.lwd = 2)
legend("topright",
       legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
       fill = c("white", rev(color.jet(11))),
       border = NA,
       y.intersp = .5,
       cex = 1, text.font = 1,
       title = "Correlation", title.adj = 0.2, title.cex = 1.25)
legend("bottomright",
       legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
       fill = node_color,
       cex = 1,
       y.intersp = .8,
       title = "Block", title.adj = 0.2, title.cex = 1.25)
text(.6, .95, labels = "I", family = "arial", cex = 2)
text(-.9, -.95, labels = "II", family = "arial", cex = 2)
dev.off() 

# S11d

matrix_il_combined_1_4 <- plot_diablo(diablo_il_combined_1_4, cutoff_value = 0.8)

g <- graph_from_adjacency_matrix(matrix_il_combined_1_4, mode = "undirected", weighted = T, diag = F)
cluster = list("1" = c("pcd_phe", "pcd_pro", "il_pro", "pcd_gly", "il_gly",
                       "p_K03074", "p_K03072", "p_K03413", "p_K07699", "p_K12257",
                       "g_K02169", "g_K05934", "g_K12984", 
                       "g_Bifidobacterium boum",
                       "ssc_396919", "ssc_100620246", "ssc_396807", "ssc_396766", 
                       "ssc_397197", "ssc_445461", "ssc_100514249", "ssc_100157995",
                       "ssc_100158244", "ssc_397316", 
                       "A0A9D4ZV85_PEA", "A0A9D4VX07_PEA", "A0A9D4XPB1_PEA",
                       "A0A9D5A4B1_PEA", "A0A9D5A4J7_PEA",
                       "m_glutamine", "m_phenylalanine", "m_valine", "m_aspartate", 
                       "m_asparagine", "m_aspartate", "m_glutamate", "m_trimethylamine", "m_acetate"),
               "2" = c("il_total_starch", "il_insp6", "il_ca", "pcd_ca", "pcd_p",
                       "p_K12686", "p_K03281", "p_K07571", "p_K02372", "p_Enterocloster",
                       "g_K00863", "g_K01961", "g_K03718",
                       "g_Hominenteromicrobium mulieris", "g_Bacteroides xylanisolvens",
                       "g_Bariatricus sp004560705", 
                       "ssc_100521982", "ssc_100624628", "ssc_397602", "ssc_396921",
                       "A0A9D4XTF6_PEA", "A0A9D5BQA3_PEA", "A0A9D4XKQ8_PEA",
                       "A0A9D5B8M5_PEA",
                       "m_trimethylamine", "m_propionate"),
               "3" = c("enz_try", "enz_chy", "enz_amy", "il_p", "il_total_starch", 
                       "p_K17810", "p_K01493", "p_K05910", "p_K03218",
                       "g_K03980", "g_K06956", "g_K02502", "g_K03718",
                       "g_Fusobacterium gastrosuis", 
                       "ssc_100155038", "ssc_808504", "ssc_397113",
                       "A0A9D5BN09_PEA", "A0A9D4YP35_PEA", "A0A9D4WQC0_PEA",
                       "A0A9D4VWB7_PEA", "A0A9D4YG61_PEA",
                       "m_lactate", "m_tyrosine"),
               "4" = c("p_K00854", "p_K01738", "p_K19168", "p_K00962",
                       "g_K00873", "g_K14261", "g_K06074",
                       "ssc_595116", 
                       "A0A9D4VQB0_PEA",
                       "m_butyrate"))
layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
                                               vcount = vcount(g), area = vcount(g)^2.95, repulse.rad=vcount(g)^3.5)
layout_fr <- norm_coords(layout_fr)
# assign colors to weight values
color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
E(g)$color <- color.jet(100)[color_index]
E(g)$width <- .75
node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                          str_detect(V(g)$name, "^g_") ~ 2,
                          str_detect(V(g)$name, "^p_") ~ 3,
                          str_detect(V(g)$name, "^ssc_") ~ 4,
                          str_detect(V(g)$name, "_PEA$") ~ 5,
                          str_detect(V(g)$name, "^m_") ~ 6,
                          .default = 0)
V(g)$color <- node_color[color_vector]
V(g)$label.family <- "arial"
V(g)$label.font <- 1
V(g)$label.color <- "black"
V(g)$label.cex <- .75
jpeg(filename = paste0("plots/", get_script_number(), "_figS11d", save_name, ".jpeg"),
     width = 20, height = 15, unit="cm", res = 500, pointsize = 12, family = "arial")
par(mar = c(0,0,0,0))
plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = F,
     vertex.shape = "rectangle", 
     vertex.size = str_width(V(g)$name)*1.8,
     vertex.size2 = 3,
     vertex.frame.width = .1,
     mark.groups = cluster,
     mark.shape = 0.5,
     mark.expand = 0.5,
     mark.col = NA,
     mark.border = "grey10",
     mark.lwd = 2)
legend("topright",
       legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
       fill = c("white", rev(color.jet(11))),
       border = NA,
       y.intersp = .5,
       cex = 1, text.font = 1,
       title = "Correlation", title.adj = 0.2, title.cex = 1.25)
legend("bottomright",
       legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
       fill = node_color,
       y.intersp = .8,
       cex = 1,
       title = "Block", title.adj = 0.2, title.cex = 1.25)
text(-.7, .2, labels = "I", family = "arial", cex = 2)
text(.6, .95, labels = "II", family = "arial", cex = 2)
text(.95, .6, labels = "III", family = "arial", cex = 2)
text(.8, -.85, labels = "IV", family = "arial", cex = 2)
dev.off() 

# S12d

matrix_fa_combined_1_2 <- plot_diablo(diablo_fa_combined_1_2, cutoff_value = 0.75)

g <- graph_from_adjacency_matrix(matrix_fa_combined_1_2, mode = "undirected", weighted = T, diag = F)
cluster = list("1" = c("enz_try", "fa_insp6", "hg_insp6",
                       "p_K02652", "p_K02243", "p_K03544", "p_K22339",
                       "p_Phascolarctobacterium", 
                       "g_K01267", "g_K02535", "g_K11732", "g_K16363", 
                       "g_Ruminococcus_C sp937915125",
                       "ssc_100521982", "ssc_406192", "ssc_733683",
                       "ssc_396685", "ssc_100312960",
                       "A0A9D4VYD1_PEA", "A0A9D5BH02_PEA", "A0A9D4Y6S0_PEA",
                       "A0A9D4VY92_PEA", "A0A9D4X2N1_PEA",
                       "m_succinate", "m_acetate"),
               "2" = c("hg_ge", "hg_cp", "fa_p",
                       "p_K00440", "p_K02032", "p_K07386", "p_K01622", "p_K04088",
                       "p_K04087", "p_K10192", "p_K03837",
                       "p_RGIG1693", "p_Hominicoprocola",
                       "g_K07068", "g_K07701", "g_K17319", "g_K01676",
                       "g_CAG-103 sp905215475",
                       "ssc_100627480", "ssc_397140", "ssc_397419",
                       "A0A9D5B2R9_PEA", "A0A9D4VRF6_PEA", "A0A068LJH6_PEA",
                       "m_3-phenylpropionate", "m_valerate", "m_propionate"),
               "3" = c("fa_ti", "fa_cp",
                       "p_K07335", "p_K07337", "p_K03832", "p_K00348", 
                       "p_PeH17", "p_UBA4363",
                       "g_K01846", "g_K01442", "g_K07016",
                       "ssc_397192", "ssc_100522855", "ssc_100154047",
                       "A0A9D4W4Y0_PEA", "A0A9D4WKA7_PEA", "D3VND9_PEA", "Q9M3X6_PEA",
                       "A0A9D4XWP6_PEA",
                       "m_valine", "m_phenylalanine", "m_methionine", "m_tyrosine"))
layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
                                               vcount = vcount(g), area = vcount(g)^3, repulse.rad=vcount(g)^3.4)
layout_fr <- norm_coords(layout_fr)
# assign colors to weight values
color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
E(g)$color <- color.jet(100)[color_index]
E(g)$width <- .75
node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                          str_detect(V(g)$name, "^g_") ~ 2,
                          str_detect(V(g)$name, "^p_") ~ 3,
                          str_detect(V(g)$name, "^ssc_") ~ 4,
                          str_detect(V(g)$name, "_PEA$") ~ 5,
                          str_detect(V(g)$name, "^m_") ~ 6,
                          .default = 0)
V(g)$color <- node_color[color_vector]
V(g)$label.family <- "arial"
V(g)$label.font <- 1
V(g)$label.color <- "black"
V(g)$label.cex <- .8
jpeg(filename = paste0("plots/", get_script_number(), "_figS12d", save_name, ".jpeg"),
     width = 20, height = 15, unit="cm", res = 500, pointsize = 12, family = "arial")
par(mar = c(0,0,0,0))
plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = F,
     vertex.shape = "rectangle", 
     vertex.size = str_width(V(g)$name)*1.9,
     vertex.size2 = 3.5,
     vertex.frame.width = .1,
     mark.groups = cluster,
     mark.shape = 0.5,
     mark.expand = 0.5,
     mark.col = NA,
     mark.border = "grey10",
     mark.lwd = 2)
legend("topright",
       legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
       fill = c("white", rev(color.jet(11))),
       border = NA,
       y.intersp = .5,
       cex = 1, text.font = 1,
       title = "Correlation", title.adj = 0.2, title.cex = 1.25)
legend("bottomright",
       legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
       fill = node_color,
       y.intersp = .8,
       cex = 1,
       title = "Block", title.adj = 0.2, title.cex = 1.25)
text(.7, .95, labels = "I", family = "arial", cex = 2)
text(1, .75, labels = "II", family = "arial", cex = 2)
text(-.9, -.9, labels = "III", family = "arial", cex = 2)
dev.off() 

# S13d

matrix_fa_combined_1_3 <- plot_diablo(diablo_fa_combined_1_3, cutoff_value = 0.8)

g <- graph_from_adjacency_matrix(matrix_fa_combined_1_3, mode = "undirected", weighted = T, diag = F)
cluster = list("1" = c("fa_p", "hg_cp", "fa_ca",
                       "p_K03231", "p_K01649", "p_K00582", "p_K00123", "p_K00757",
                       "p_K01596", "p_K00609", "p_K14126", "p_K21990", "p_K01191",
                       "p_Clostridium_AI", 
                       "g_K03420", "g_K14098", "g_K07558", "g_K16792", "g_K14115",
                       "ssc_397108", "ssc_397376", "ssc_100521982", 
                       "A0A9D5B9N7_PEA", 
                       "m_acetate", "m_valerate", "m_uracil", "m_propionate", 
                       "m_butyrate", "m_aspartate", "m_isobutyrate"),
               "2" = c("fa_cp",
                       "p_K18676", "p_K01915", "p_K02035", "p_K02051",
                       "p_PeH17",
                       "g_K01031", "g_K04835", "g_K01846", "g_K00772", "g_K01758",
                       "ssc_396921", "ssc_445518", "ssc_396674", "ssc_445532",
                       "A0A9D4YG61_PEA", "A0A9D4W4Y0_PEA", "A0A9D4W9W0_PEA",
                       "A0A9D5AY78_PEA", "A0A9D5A4B1_PEA",
                       "m_methionine", "m_phenylalanine", "m_tyrosine", "m_valine"))
layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
                                               vcount = vcount(g), area = vcount(g)^3, repulse.rad=vcount(g)^3.5)
layout_fr <- norm_coords(layout_fr)
# assign colors to weight values
color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
E(g)$color <- color.jet(100)[color_index]
E(g)$width <- .75
node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                          str_detect(V(g)$name, "^g_") ~ 2,
                          str_detect(V(g)$name, "^p_") ~ 3,
                          str_detect(V(g)$name, "^ssc_") ~ 4,
                          str_detect(V(g)$name, "_PEA$") ~ 5,
                          str_detect(V(g)$name, "^m_") ~ 6,
                          .default = 0)
V(g)$color <- node_color[color_vector]
V(g)$label.family <- "arial"
V(g)$label.font <- 1
V(g)$label.color <- "black"
V(g)$label.cex <- .9
jpeg(filename = paste0("plots/", get_script_number(), "_figS13d", save_name, ".jpeg"),
     width = 20, height = 15, unit="cm", res = 500, pointsize = 12, family = "arial")
par(mar = c(0,0,0,0))
plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = F,
     vertex.shape = "rectangle", 
     vertex.size = str_width(V(g)$name)*2.1,
     vertex.size2 = 3.5,
     vertex.frame.width = .1,
     mark.groups = cluster,
     mark.shape = 0.5,
     mark.expand = 0.5,
     mark.col = NA,
     mark.border = "grey10",
     mark.lwd = 2)
legend("topright",
       legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
       fill = c("white", rev(color.jet(11))),
       border = NA,
       y.intersp = .5,
       cex = 1, text.font = 1,
       title = "Correlation", title.adj = 0.2, title.cex = 1.25)
legend("bottomright",
       legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
       fill = node_color,
       y.intersp = .8,
       cex = 1,
       title = "Block", title.adj = 0.2, title.cex = 1.25)
text(.7, .95, labels = "I", family = "arial", cex = 2)
text(-.9, -.9, labels = "II", family = "arial", cex = 2)
dev.off() 


# redefine functions for kegg-taxa origin 

proteins_func_tax_combined <- readRDS("clean/4_function_taxa_combined_long_norm_imp_rel_abd_filtered.RDS") %>%
  inner_join(meta, by = "sampleid")
proteins_func_tax_combined_il <- filter_ileum(proteins_func_tax_combined)
proteins_func_tax_combined_fa <- filter_faeces(proteins_func_tax_combined)

p_plot_kegg_taxa_origin <- function(input_df, kegg_ko_filter, taxa_level = "G", threshold = 0.01) {
  filtered_df <- input_df %>%
    filter(str_detect(kegg_ko, kegg_ko_filter)) %>%
    group_by(!!sym(taxa_level), proteinid) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
    group_by(!!sym(taxa_level)) %>%
    summarise(rel_abd = sum(rel_abd)) %>%
    mutate(rel_rel_abd = rel_abd / sum(rel_abd)) %>%
    ungroup() %>%
    mutate(!!sym(taxa_level) := ifelse(rel_rel_abd < threshold, "other", !!sym(taxa_level))) %>%
    group_by(!!sym(taxa_level)) %>%
    summarize(rel_abd = sum(rel_abd), .groups = "drop") 
  filtered_df_colored <- filtered_df %>%
    left_join(color_df, by = c("G" = "name")) %>%
    add_column("color2" = sample(colors, size = nrow(filtered_df))) %>%
    mutate(color = ifelse(is.na(color), color2, color))

  if ("other" %in% filtered_df[[taxa_level]]) {
    filtered_df <- filtered_df %>%
      mutate(!!sym(taxa_level) := factor(!!sym(taxa_level), levels = c(unique(!!sym(taxa_level))[-which(unique(!!sym(taxa_level))=="other")], "other")))
  } else {
    filtered_df <- filtered_df %>%
      mutate(!!sym(taxa_level) := factor(!!sym(taxa_level), levels = !!sym(taxa_level)))
  }

  p <- ggplot(filtered_df, aes(x = 1, y = rel_abd, fill = !!sym(taxa_level))) +
    geom_bar(stat = "identity", position = "stack", width = .6) +
    scale_fill_manual(values = filtered_df_colored$color) +
    scale_y_continuous(limits = c(0, sum(filtered_df$rel_abd)*1.05), expand = c(0,0)) +
    labs(x = kegg_ko_filter, y = "Relative abundance (%)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(p)
}

p_plot_kegg_for_taxa <- function(input_df, taxa_name, threshold = 1, use_names = F) {
  filtered_table <- input_df %>%
    filter(G == taxa_name) %>%
    calculate_kegg_abundance(abundance_column = "rel_abd") %>%
    group_by(sampleid) %>%
    mutate(rel_abd = rel_abd / sum(rel_abd) * 100) %>%
    ungroup() %>%
    left_join(meta, by = "sampleid") %>%
    group_by(kegg_ko) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
    mutate(kegg_ko = ifelse(rel_abd < threshold, "other", kegg_ko)) %>%
    group_by(kegg_ko) %>%
    summarise(rel_abd = sum(rel_abd), .groups = "drop") %>%
    mutate(kegg_ko = factor(kegg_ko, levels = c(unique(kegg_ko)[-which(unique(kegg_ko)=="other")], "other")))
  
  if (isTRUE(use_names)) {
    filtered_table_ko <- annotate_keggs(unique(filtered_table$kegg_ko))
    
    filtered_table_kegg <- left_join(filtered_table, filtered_table_ko, by = "kegg_ko")
    
    p <- ggplot(filtered_table_kegg, aes(x = 1, y = rel_abd, fill = reorder(name, rel_abd))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      scale_y_continuous(limits = c(0, 100.01), expand = c(0,0)) +
      labs(x = taxa_name, fill = "KEGG ortholog", y = "Relative abundance (%)", x = "") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- ggplot(filtered_table, aes(x = 1, y = rel_abd, fill = reorder(kegg_ko, rel_abd))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      scale_y_continuous(limits = c(0, 100.01), expand = c(0,0)) +
      labs(x = taxa_name, fill = "KEGG ortholog", y = "Relative abundance (%)", x = "") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  return(p)
}

# figure 4
# to be determined

# matrix_il_combined_1_2 <- plot_diablo(diablo_il_combined_1_2, cutoff_value = 0.8)
# 
# pdf(NULL)
# matrix_il_combined_1_2 <- circosPlot(diablo_il_combined_1_2, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
# dev.off()
# 
# matrix_il_combined_1_2_sub1sub <- filter_submatrix(matrix_il_combined_1_2,
#                                          vector = c( "il_tyr", "il_ile", "il_ala",
#                                                     "il_val", "il_ser", "il_leu", "il_his", "il_asp", "il_lys",
#                                                     "p_K05910", "p_K00700", 
#                                                     "g_Prevotella pectinovora",
#                                                     "g_Prevotella copri",
#                                                     "g_Prevotella copri_A",
#                                                     "ssc_733607"))
# 
# g <- graph_from_adjacency_matrix(matrix_il_combined_1_2_sub1sub, mode = "undirected", weighted = T, diag = F)
# layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
# layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
#                                                vcount = vcount(g), area = vcount(g)^3, repulse.rad=vcount(g)^3.4)
# layout_fr <- norm_coords(layout_fr)
# # assign colors to weight values
# color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
# E(g)$color <- color.jet(100)[color_index]
# E(g)$width <- 1.75
# node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
# color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
#                           str_detect(V(g)$name, "^g_") ~ 2,
#                           str_detect(V(g)$name, "^p_") ~ 3,
#                           str_detect(V(g)$name, "^ssc_") ~ 4,
#                           str_detect(V(g)$name, "_PEA$") ~ 5,
#                           str_detect(V(g)$name, "^m_") ~ 6,
#                           .default = 0)
# V(g)$color <- node_color[color_vector]
# V(g)$label.family <- "arial"
# V(g)$label.font <- 1
# V(g)$label.color <- "black"
# V(g)$label.cex <- 1
# jpeg(filename = paste0("plots/", get_script_number(), "_figS4a", save_name, ".jpeg"),
#      width = 10, height = 7.5, unit="cm", res = 500, pointsize = 12, family = "arial")
# par(mar = c(0,0,0,0))
# plot(g, layout = layout.circle, 
#      margin = c(0,0,0,0), rescale = F,
#      vertex.shape = "rectangle", 
#      vertex.size = str_width(V(g)$name)*4.5,
#      vertex.size2 = 7.5,
#      vertex.frame.width = .1)
# # legend("topright",
# #        legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
# #        fill = c("white", rev(color.jet(11))),
# #        border = NA,
# #        y.intersp = .5,
# #        cex = 1, text.font = 1,
# #        title = "Correlation", title.adj = 0.2, title.cex = 1.25)
# # legend("bottomright",
# #        legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
# #        fill = node_color,
# #        y.intersp = .8,
# #        cex = 1,
# #        title = "Block", title.adj = 0.2, title.cex = 1.25)
# dev.off() 
# fig4a <- as.ggplot(pheatmap(matrix_il_combined_1_2_sub1sub))

matrix_il_combined_1_2 <- plot_diablo(diablo_il_combined_1_2, cutoff_value = 0.875)

g <- graph_from_adjacency_matrix(matrix_il_combined_1_2, mode = "undirected", weighted = T, diag = F)
cluster <- list("1" = c("il_met", "il_tyr", "il_phe", "il_thr", "il_ile", "il_ala",
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
                        "m_tyrosine", "m_valine", "m_phenylalanine"),
                "2" = c("il_tdf",
                        "p_K01200", "p_K01582", "p_K01585", "p_K01581",
                        "g_K02444", "g_K11782", "g_K02355", "g_K09949", 
                        "g_Streptococcus oriscaviae", "g_Mitsuokella multacida",
                        "ssc_808504", "ssc_733685", "ssc_397113",
                        "A0A9D5BQA3_PEA",
                        "m_propionate", "m_trimethylamine"),
                "3" = c("il_ip_4_1256", 
                        "p_K11929", "p_K09475", "p_K09476", "p_K14062", "p_K16076", "p_K01533",
                        "g_K01496", "g_K13935",
                        "g_Nanosynbacter sp029975625", "g_Escherichia sp004211955",
                        "g_Flavobacterium psychrophilum_B", "g_Escherichia sp005843885",
                        "g_JALHET01 sp022839525", "g_Advenella sp023423975",
                        "A0A9D4WQ37_PEA"), 
                "4" = c("m_galactose",
                        "p_K15582", "p_K02913", "p_K03706", "p_K02034", "p_K03768",
                        "ssc_100158117"),
                "5" = c("pcd_ca",
                        "p_K00020", "p_K00042", "p_K03336", 
                        "g_K00012", "g_K06871",
                        "ssc_654405", "ssc_100152549",
                        "A0A9D4WM15_PEA"))
layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
                                               vcount = vcount(g), area = vcount(g)^3, repulse.rad=vcount(g)^3.65)
layout_fr <- norm_coords(layout_fr)
# assign colors to weight values
color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
E(g)$color <- color.jet(100)[color_index]
E(g)$width <- .75
node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                          str_detect(V(g)$name, "^g_") ~ 2,
                          str_detect(V(g)$name, "^p_") ~ 3,
                          str_detect(V(g)$name, "^ssc_") ~ 4,
                          str_detect(V(g)$name, "_PEA$") ~ 5,
                          str_detect(V(g)$name, "^m_") ~ 6,
                          .default = 0)
V(g)$color <- node_color[color_vector]
V(g)$label.family <- "arial"
V(g)$label.font <- 1
V(g)$label.color <- "black"
V(g)$label.cex <- .75
jpeg(filename = paste0("plots/", get_script_number(), "_fig4a", save_name, ".jpeg"),
     width = 20, height = 15, unit="cm", res = 500, pointsize = 12, family = "arial")
par(mar = c(0,0,0,0))
plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = F,
     vertex.shape = "rectangle", 
     vertex.size = str_width(V(g)$name)*1.8,
     vertex.size2 = 3,
     vertex.frame.width = .1, 
     mark.groups = cluster,
     mark.shape = .5,
     mark.expand = 0.5, 
     mark.col = NA,
     mark.border = "grey10",
     mark.lwd = 2)
legend("topright",
       legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
       fill = c("white", rev(color.jet(11))),
       border = NA,
       y.intersp = .5,
       cex = 1, text.font = 1,
       title = "Correlation", title.adj = 0.2, title.cex = 1.25)
legend("bottomright",
       legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
       fill = node_color,
       y.intersp = .8,
       cex = 1,
       title = "Block", title.adj = 0.2, title.cex = 1.25)
text(.9, .95, labels = "I", family = "arial", cex = 2)
text(.9, -.95, labels = "II", family = "arial", cex = 2)
text(-.7, -.95, labels = "III", family = "arial", cex = 2)
text(-1, -.5, labels = "IV", family = "arial", cex = 2)
text(-1, .5, labels = "V", family = "arial", cex = 2)
dev.off() 

pdf(NULL)
matrix_il_combined_1_2 <- circosPlot(diablo_il_combined_1_2, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

fig4b <- as.ggplot(pheatmap(il_combined_1_2_sub1, fontsize_row = 8, fontsize_col = 14))

empty <- ggplot() +
  theme_void()

fig4 <- (empty / fig4b) +
  plot_layout(heights = c(4,4)) +
  plot_annotation(tag_levels = list(c("a", "b"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_fig4.jpeg",
       plot = fig4,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300)

# Fig 5

matrix_fa_combined_1_4 <- plot_diablo(diablo_fa_combined_1_4, cutoff_value = 0.8)

g <- graph_from_adjacency_matrix(matrix_fa_combined_1_4, mode = "undirected", weighted = T, diag = F)
cluster = list("1" = c("fa_k", "hg_dm", "hg_tdf", "fa_ti", 
                        "p_K06410", "p_K04079", "p_K02652", "p_K02662", "p_K02243",
                        "g_K01267",
                        "g_UMGS124 sp019420325", "g_UMGS124 sp902464015", 
                        "g_HGM13006 sp029012465", "g_UMGS124 sp900555105",
                        "ssc_100525899", "ssc_407610", "ssc_100037943", "ssc_445461",
                        "ssc_397397",
                        "A0A9D5BQA3_PEA", "A0A9D5BN09_PEA", "A0A9D4XTC5_PEA", 
                        "A0A9D4VJF2_PEA", "A0A9D4WCI6_PEA",
                        "m_acetate", "m_butyrate", "m_propionate", "m_methionine",
                        "m_galactose"),
                "2" = c("enz_carb", "enz_try", "enz_chy", 
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
layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
layout_fr <- qgraph.layout.fruchtermanreingold(get.edgelist(g, names = F), #weights = 1-(abs(E(g)$weight)), 
                                               vcount = vcount(g), area = vcount(g)^3, repulse.rad=vcount(g)^3.54)
layout_fr <- norm_coords(layout_fr)
# assign colors to weight values
color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
E(g)$color <- color.jet(100)[color_index]
E(g)$width <- .75
node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                          str_detect(V(g)$name, "^g_") ~ 2,
                          str_detect(V(g)$name, "^p_") ~ 3,
                          str_detect(V(g)$name, "^ssc_") ~ 4,
                          str_detect(V(g)$name, "_PEA$") ~ 5,
                          str_detect(V(g)$name, "^m_") ~ 6,
                          .default = 0)
V(g)$color <- node_color[color_vector]
V(g)$label.family <- "arial"
V(g)$label.font <- 1
V(g)$label.color <- "black"
V(g)$label.cex <- .8
jpeg(filename = paste0("plots/", get_script_number(), "_fig5a", save_name, ".jpeg"),
     width = 20, height = 15, unit="cm", res = 500, pointsize = 12, family = "arial")
par(mar = c(0,0,0,0))
plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = F,
     vertex.shape = "rectangle", 
     vertex.size = str_width(V(g)$name)*2.0,
     vertex.size2 = 3.5,
     vertex.frame.width = .1,
     mark.groups = cluster,
     mark.shape = 0.5,
     mark.expand = 0.5,
     mark.col = NA,
     mark.border = "grey10",
     mark.lwd = 2)
legend("topright",
       legend = c(NA,1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
       fill = c("white", rev(color.jet(11))),
       border = NA,
       y.intersp = .5,
       cex = 1, text.font = 1,
       title = "Correlation", title.adj = 0.2, title.cex = 1.25)
legend("bottomright",
       legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
       fill = node_color,
       y.intersp = .8,
       cex = 1,
       title = "Block", title.adj = 0.2, title.cex = 1.25)
text(.6, .95, labels = "I", family = "arial", cex = 2)
text(-.9, -.95, labels = "II", family = "arial", cex = 2)
dev.off() 

pdf(NULL)
matrix_fa_combined_1_4 <- circosPlot(diablo_fa_combined_1_4, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

fig5b <- as.ggplot(pheatmap(fa_combined_1_4_sub1, fontsize_row = 12, fontsize_col = 16))

fig5c <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_fa, kegg_ko_filter = "K02662") +
  labs(fill = "Genus", title = "K02662") +
  theme(axis.title.x = element_blank())

empty <- ggplot() +
  theme_void()

fig5bc <- (fig5b | fig5c) +
  plot_layout(widths = c(24,2))

fig5 <- (empty / fig5bc) +
  plot_layout(heights = c(4,4)) +
  plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_fig5.jpeg",
       plot = fig5,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300) 

# Heat maps
library(ggplotify)

pdf(NULL)
matrix_il_combined_1_2 <- circosPlot(diablo_il_combined_1_2, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

il_combined_1_2_sub2 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("il_tdf",
                                                    "p_K01200", "p_K01582", "p_K01585", "p_K01581",
                                                    "g_K02444", "g_K11782", "g_K02355", "g_K09949", 
                                                    "g_Streptococcus oriscaviae", "g_Mitsuokella multacida",
                                                    "ssc_808504", "ssc_733685", "ssc_397113",
                                                    "A0A9D5BQA3_PEA",
                                                    "m_propionate", "m_trimethylamine"))

il_combined_1_2_sub3 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("il_ip_4_1256", 
                                                    "p_K11929", "p_K09475", "p_K09476", "p_K14062", "p_K16076",
                                                    "g_K01496",
                                                    "g_Nanosynbacter sp029975625", "g_Escherichia sp004211955",
                                                    "g_Flavobacterium psychrophilum_B", "g_Escherichia sp005843885",
                                                    "g_JALHET01 sp022839525", "g_Advenella sp023423975",
                                                    "A0A9D4WQ37_PEA"))

il_combined_1_2_sub4 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("m_galactose",
                                                    "p_K15582", "p_K02913", "p_K03706", "p_K02034",
                                                    "ssc_100158117"))

il_combined_1_2_sub5 <- filter_submatrix(matrix_il_combined_1_2,
                                         vector = c("pcd_ca",
                                                    "p_K00020", "p_K00042", "p_K03336", 
                                                    "g_K00012", "g_K06871",
                                                    "ssc_654405", "ssc_100152549",
                                                    "A0A9D4WM15_PEA"))



#a <- as.ggplot(pheatmap(il_combined_1_2_sub1, fontsize_row = 10))
a <- as.ggplot(pheatmap(il_combined_1_2_sub2, fontsize = 16))
# c <- as.ggplot(pheatmap(il_combined_1_2_sub3))
# d <- as.ggplot(pheatmap(il_combined_1_2_sub4))
# e <- as.ggplot(pheatmap(il_combined_1_2_sub5))

b <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_il, kegg_ko_filter = "K00700") +
  labs(fill = "Genus")

c <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_il, kegg_ko_filter = "K01200") +
  labs(fill = "Genus")

figS15 <- (a | (b/c)) +
  plot_layout(widths = c(12,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))
  
ggsave(filename= "92_figS15.jpeg",
       plot = figS15,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 12,
       scale=2,
       dpi=300)         


pdf(NULL)
matrix_il_combined_1_3 <- circosPlot(diablo_il_combined_1_3, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

il_combined_1_3_sub2 <- filter_submatrix(matrix_il_combined_1_3,
                                         vector = c("il_total_starch",
                                                    "p_K06305", "p_K03640", "p_K01176", "p_K02217",
                                                    "p_Basfia_A",
                                                    "g_K10542", "g_K03775", "g_K14441", "g_K03592", "g_K04085",
                                                    "ssc_397060", "ssc_100511536", "ssc_808504", "ssc_397113",
                                                    "A0A9D4VWB7_PEA", "A0A9D4VWB7_PEA", "A0A9D4ZWX1_PEA",
                                                    "A0A9D5BJS9_PEA", "A0A9D4WDI6_PEA",
                                                    "m_trimethylamine"))

a <- as.ggplot(pheatmap(il_combined_1_3_sub1))
b <- as.ggplot(pheatmap(il_combined_1_3_sub2, fontsize_col = 14, fontsize_row =10))

c <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_il, kegg_ko_filter = "K07699") +
  labs(fill = "Genus")

d <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_il, kegg_ko_filter = "K01176") +
  labs(fill = "Genus")

cd <- (c / d) +
  plot_layout(heights = c(4,4))

bc <- (b | cd) +
  plot_layout(widths = c(12,4))

figS16 <- a /
  bc +
  plot_layout(heights = c(8,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS16.jpeg",
       plot = figS16,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300) 

pdf(NULL)
matrix_il_combined_1_4 <- circosPlot(diablo_il_combined_1_4, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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
                                                    "m_asparagine", "m_aspartate", "m_glutamate", "m_trimethylamine",
                                                    "m_acetate"))

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

il_combined_1_4_sub3 <- filter_submatrix(matrix_il_combined_1_4,
                                         vector = c("enz_try", "enz_chy", "enz_amy", "il_p", "il_total_starch", 
                                                    "p_K17810", "p_K01493", "p_K05910", "p_K03218",
                                                    "g_K03980", "g_K06956", "g_K02502", "g_K03718",
                                                    "g_Fusobacterium gastrosuis", 
                                                    "ssc_100155038", "ssc_808504", "ssc_397113",
                                                    "A0A9D5BN09_PEA", "A0A9D4YP35_PEA", "A0A9D4WQC0_PEA",
                                                    "A0A9D4VWB7_PEA", "A0A9D4YG61_PEA",
                                                    "m_lactate", "m_tyrosine"))

il_combined_1_4_sub4 <- filter_submatrix(matrix_il_combined_1_4,
                                         vector = c("p_K00854", "p_K01738", "p_K19168", "p_K00962",
                                                    "g_K00873", "g_K14261", "g_K06074",
                                                    "ssc_595116", 
                                                    "A0A9D4VQB0_PEA",
                                                    "m_butyrate"))

a <- as.ggplot(pheatmap(il_combined_1_4_sub1, fontsize = 16))
b <- as.ggplot(pheatmap(il_combined_1_4_sub2))
c <- as.ggplot(pheatmap(il_combined_1_4_sub3))
d <- as.ggplot(pheatmap(il_combined_1_4_sub4))

c <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_il, kegg_ko_filter = "K01176") +
  labs(fill = "Genus")


figS17 <- (a | b) /
  (c | d) +
  plot_layout(heights = c(5,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS17.jpeg",
       plot = a,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 20,
       scale=2,
       dpi=300) 

pdf(NULL)
matrix_fa_combined_1_2 <- circosPlot(diablo_fa_combined_1_2, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

fa_combined_1_2_sub2 <- filter_submatrix(matrix_fa_combined_1_2,
                                         vector = c("fa_ti", "fa_cp",
                                                    "p_K07335", "p_K07337", "p_K03832", "p_K00348", 
                                                    "p_PeH17", "p_UBA4363",
                                                    "g_K01846", "g_K01442", "g_K07016",
                                                    "ssc_397192", "ssc_100522855", "ssc_100154047",
                                                    "A0A9D4W4Y0_PEA", "A0A9D4WKA7_PEA", "D3VND9_PEA", "Q9M3X6_PEA",
                                                    "A0A9D4XWP6_PEA",
                                                    "m_valine", "m_phenylalanine", "m_methionine", "m_tyrosine"))

a <- as.ggplot(pheatmap(fa_combined_1_2_sub3, fontsize_row = 14))
c <- as.ggplot(pheatmap(fa_combined_1_2_sub2, fontsize = 14))
#b <- as.ggplot(pheatmap(fa_combined_1_2_sub1))

b <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_fa, kegg_ko_filter = "K22339") +
  labs(fill = "Genus")

d <- p_plot_kegg_for_taxa(input_df = proteins_func_tax_combined_fa, taxa_name = "PeH17", use_names = F)


ab <- (a | b) +
  plot_layout(widths = c(18,4))

cd <- (c | d) +
  plot_layout(widths = c(18,4))


figS18 <- (ab / cd) +
  plot_layout(heights = c(4,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS18.jpeg",
       plot = figS18,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300) 

pdf(NULL)
matrix_fa_combined_1_3 <- circosPlot(diablo_fa_combined_1_3, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

fa_combined_1_3_sub2 <- filter_submatrix(matrix_fa_combined_1_3,
                                         vector = c("fa_cp",
                                                    "p_K18676", "p_K01915", "p_K02035", "p_K02051",
                                                    "p_PeH17",
                                                    "g_K01031", "g_K04835", "g_K01846", "g_K00772", "g_K01758",
                                                    "ssc_396921", "ssc_445518", "ssc_396674", "ssc_445532",
                                                    "A0A9D4YG61_PEA", "A0A9D4W4Y0_PEA", "A0A9D4W9W0_PEA",
                                                    "A0A9D5AY78_PEA", "A0A9D5A4B1_PEA",
                                                    "m_methionine", "m_phenylalanine", "m_tyrosine", "m_valine"))

a <- as.ggplot(pheatmap(fa_combined_1_3_sub1, fontsize_row = 14))
c <- as.ggplot(pheatmap(fa_combined_1_3_sub2, fontsize = 14))

b <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_fa, kegg_ko_filter = "K01596") +
  labs(fill = "Genus")

ab <- (a | b) +
  plot_layout(widths = c(16,4))

figS19 <- ab / c +
  plot_layout(heights = c(4,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS19.jpeg",
       plot = figS19,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 25,
       scale=2,
       dpi=300) 

pdf(NULL)
matrix_fa_combined_1_4 <- circosPlot(diablo_fa_combined_1_4, cutoff = 0, size.variables = 1, line = T, size.labels = 1.5)
dev.off()

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

fa_combined_1_4_sub2 <- filter_submatrix(matrix_fa_combined_1_4,
                                         vector = c("enz_carb", "enz_try", "enz_chy",
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

#a <- as.ggplot(pheatmap(fa_combined_1_4_sub1, fontsize = 12))
a <- as.ggplot(pheatmap(fa_combined_1_4_sub2, fontsize_row = 12))

b <- p_plot_kegg_taxa_origin(input_df = proteins_func_tax_combined_fa, kegg_ko_filter = "K02662") +
  labs(fill = "Genus")

figS20 <- (a | b) +
  plot_layout(widths = c(24,4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(filename= "92_figS20.jpeg",
       plot = figS20,
       device= "jpeg", 
       path = "plots", 
       units = "cm", 
       width = 20,
       height = 12,
       scale=2,
       dpi=300) 
