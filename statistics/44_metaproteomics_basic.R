library(here)
library(venn)

source(here("30_omics_functions.R"))
source("E:/R/source/ggplot2_theme_bw.R")

# load meta

meta <- readRDS("clean/meta1.RDS")

# load summary file

summary_long <- readRDS("clean/4_summary_long.RDS") %>%
  filter(!parameter == "ms/ms") %>%
  left_join(dplyr::select(meta, sampleid, matrix, animal, period, diet), by = "sampleid") %>%
  mutate(parameter = case_when(parameter == "ms/ms_identified" ~ "MS identified",
                               parameter == "ms/ms_identified_[%]" ~ "MS identified %",
                               parameter == "peptide_sequences_identified" ~ "Peptides"))

# load raw protein file 

proteins <- readRDS("clean/4_proteins_long_raw_intensity.RDS")

proteins_rel_abd <- readRDS("clean/4_proteins_long_raw_rel_abd.RDS")

# functions file for joining

functions <- readRDS("clean/4_functions_long_raw.RDS") %>%
  mutate(cog_cat_name = str_c(cog_category, cog_name, sep = " - "))

# plot

ggplot(summary_long, aes(x = sampleid, y = value, fill = matrix)) +
  geom_bar(stat = "identity", show.legend = F) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(rows = vars(parameter), scales = "free_y") +
  labs(y = "") +
  scale_fill_manual(values = c("grey50", "grey70")) + 
  theme(axis.text.x = element_text(size = 8))

save_big("44_ms_summary_bar")

# plot as boxplots

ggplot(summary_long, aes(x = factor(matrix, levels = c("ileal digesta", "faeces")), y = value)) +
  geom_violin(draw_quantiles = c(0.5), fill = "snow2") +
  geom_quasirandom(aes(shape = diet, color = animal), width = 0.3, size = 2) +
  #geom_jitter(aes(shape = period, color = diet), width = 0.3, size = 2) +
  facet_wrap(~parameter, scales = "free_y") +
  scale_color_manual(values = colors) +
  labs(x = "", y = "") + 
  scale_x_discrete(guide = guide_axis(angle = 30))
  

save_big("44_ms_summary_violin")

# plot number of proteingroups per sample

proteingroups_per_sample <- proteins %>%
  mutate(abu = ifelse(intensity > 0, 1, 0)) %>%
  group_by(sampleid) %>%
  summarize(abu = sum(abu)) %>%
  left_join(meta, by = "sampleid")

ggplot(proteingroups_per_sample, aes(x = factor(matrix, levels = c("ileal digesta", "faeces")), y = abu)) +
  geom_boxplot() +
  geom_jitter(width = 0.3, size = 2) +
  labs(x = "matrix", y = "number of protein groups", title = "Number of protein groups quantified")

save_big("44_proteingroups_boxplot")

# venn diagrams of protein groups

proteins_venn <- proteins %>%
  left_join(dplyr::select(meta, sampleid, matrix), by = "sampleid") %>%
  group_by(matrix, proteinid) %>%
  summarize(intensity = sum(intensity), .groups = "drop") %>%
  mutate(abu = ifelse(intensity > 0, 1, 0)) %>%
  dplyr::select(proteinid, matrix, abu) %>%
  pivot_wider(names_from = matrix, values_from = abu) %>%
  mutate(faeces = ifelse(is.na(faeces), 0, faeces),
         `ileal digesta` = ifelse(is.na(`ileal digesta`), 0, `ileal digesta`)) %>% # account for missing proteinids in different matrices
  dplyr::select(-proteinid) 

jpeg(filename = paste0("plots/", "44_proteingroups_venn", ".jpeg"), width = 7.7, height = 7.7, unit="cm", res = 1000)
venn(proteins_venn, 
     zcolor = colors,
     ilabels = "counts",
     ilcs = 0.6,
     sncs = 0.8,
     box = FALSE)
dev.off()

# plot host - micro ratio

proteins_origin <- proteins %>%
  filter(origin %in% c("pig", "micro")) %>% # only host and microbial proteins
  group_by(sampleid, origin) %>%
  summarize(intensity = sum(intensity), .groups = "drop") %>%
  group_by(sampleid) %>%
  mutate(rel_abd = intensity/sum(intensity)) %>%
  ungroup() %>%
  left_join(meta, by = "sampleid")


proteins_origin %>%
  group_by(matrix, origin) %>%
  summarize(rel_abd = mean(rel_abd), .groups = "drop") %>%
  ggplot(aes(x = factor(matrix, levels = c("ileal digesta", "faeces")), y = rel_abd, fill = origin)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_boxplot(data = filter(proteins_origin, origin == "pig"), 
               aes(x = factor(matrix, levels = c("ileal digesta", "faeces")), y = rel_abd), 
               alpha = 0, width = 0.5, lwd = 1,
               inherit.aes = F) + 
  geom_quasirandom(data = filter(proteins_origin, origin == "pig"), 
              aes(x = factor(matrix, levels = c("ileal digesta", "faeces")), y = rel_abd),
              width = 0.2, alpha = 0.5, show.legend = F) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "relative abundance", title = "Proportion of microbial and host proteins")

save_big("44_micro_host_ratio_bar")

proteins_origin %>%
  filter_ileum() %>%
  filter(origin == "micro") %>%
  arrange(desc(rel_abd)) # animal 8 often high percentage of microbial proteins

# host-micro ratio

proteins_origin_ratio <- proteins_origin %>%
  dplyr::select(-intensity) %>%
  pivot_wider(values_from = "rel_abd", names_from = "origin") %>%
  mutate(ratio = micro/pig)

combined_comparison(select_response(filter(proteins_origin_ratio, matrix == "ileal digesta"), "ratio"), 
                    transformation = "test") # no sig difference in ileum

combined_comparison(select_response(filter(proteins_origin_ratio, matrix == "faeces"), "ratio"), 
                    transformation = "test") # sig lower ratio in diet4

ggplot(filter(proteins_origin_ratio, matrix == "faeces"), aes(x = diet, y = ratio)) +
  geom_boxplot() +
  geom_quasirandom()

save_big("44_micro_host_ratio_faeces")

