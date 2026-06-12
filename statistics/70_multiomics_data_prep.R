library(here)
library(mixOmics)
library(SmCCNet)
library(DESeq2)
library(compositions)

source(here("60_correlation_functions.R"))
source(here("30_omics_functions.R"))
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
  dplyr::select(sampleno, diet, animal, period)

meta_nutrition <- meta %>%
  dplyr::select(animal, period, sampleno) %>%
  distinct()

# nutrition

pcd <- read_in_nutrition(here("data/nutrition_pcd.txt")) %>%
  prepare_nutrition_for_correlation() %>%
  dplyr::select(-c(insp6_p, n, tdf)) %>%
  dplyr::rename(tdf = tdf_single) %>%
  dplyr::rename_with(~paste0("pcd_",.x), !starts_with("sampleno"))

hindgut <- read_in_nutrition(here("data/nutrition_hindgut.txt"))  %>%
  prepare_nutrition_for_correlation() %>%
  dplyr::select(-c(insp6_p, n)) %>%
  dplyr::rename_with(~paste0("hg_",.x), !starts_with("sampleno"))

digesta_analysis <- read_in_nutrition(here("data/nutrition_digesta_analysis.txt")) %>%
  prepare_nutrition_for_correlation() %>%
  dplyr::select(-c(ip_3_126_145_245, ip_4_1234)) %>%
  dplyr::select(-c(insp6_p, n, ti...12, ti...31)) %>%
  dplyr::rename_with(~paste0("il_",.x), !starts_with("sampleno"))

faeces_analysis <- read_in_nutrition(here("data/nutrition_faeces_analysis.txt")) %>%
  prepare_nutrition_for_correlation() %>%
  dplyr::select(-total_starch) %>%
  dplyr::select(-starts_with("ip")) %>%
  dplyr::select(-c(insp6_p, n)) %>%
  dplyr::rename_with(~paste0("fa_",.x), !starts_with("sampleno"))

# enzymes

enzymes <- readRDS("clean/16_enzymes_dm.RDS")

enzymes_il <- meta %>%
  inner_join(enzymes, by = "sampleid") %>%
  filter(matrix == "ileal digesta") %>%
  dplyr::select(-c(sampleid, matrix, sampling_time, animal, period, square, diet, description, ph, dry_matter),
                -ends_with("fm")) %>%
  arrange(sampleno) %>% 
  dplyr::rename_with(~str_remove(.x, "_dm$"), !starts_with("sampleno")) %>%
  dplyr::rename_with(~paste0("enz_",.x), !starts_with("sampleno"))

enzymes_fa <- meta %>%
  inner_join(enzymes, by = "sampleid") %>%
  filter(matrix == "faeces") %>%
  dplyr::select(-c(sampleid, matrix, sampling_time, animal, period, square, diet, description, ph, dry_matter),
                -ends_with("fm")) %>%
  arrange(sampleno) %>% 
  dplyr::rename_with(~str_remove(.x, "_dm$"), !starts_with("sampleno")) %>%
  dplyr::rename_with(~paste0("enz_",.x), !starts_with("sampleno"))

# load metagenomics data as reads

breport_reads_filtered <- readRDS("clean/3_breport_reads_long_reads_filtered.RDS") %>% filter(sampleid != "103")

# kegg_reads_filtered <- readRDS("clean/3_kegg_bins_long_reads_filtered.RDS") %>% filter(sampleid != "103") # replace by contigs
kegg_reads_filtered <- readRDS("clean/3_kegg_contigs_long_reads_filtered.RDS") %>% filter(sampleid != "103")

# load log2 files of metaproteomics

proteins_norm_imp_log2 <- readRDS("clean/4_proteins_long_norm_imp_log2.RDS") %>% filter(sampleid != "103")

host_kegg_norm_imp_log2 <- readRDS("clean/4_host_kegg_long_norm_imp_log2.RDS") %>% filter(sampleid != "103")

kegg_norm_imp_log2 <- readRDS("clean/4_kegg_long_norm_imp_log2.RDS") %>% filter(sampleid != "103")

taxonomy_norm_imp_intensity <- readRDS("clean/4_taxonomy_long_norm_imp_intensity.RDS") %>% filter(sampleid != "103")

# load metabolomics data

nmr_il <- readRDS("clean/5_nmr_ileum_wide_dm.RDS") %>%
  filter(sampleid != "103") %>%
  dplyr::rename(sampleno = sampleid) %>%
  mutate(sampleno = str_extract(sampleno, "[0-9]{2}$")) %>% #prepare for correlation
dplyr::rename_with(~paste0("m_",.x), !starts_with("sampleno"))

nmr_fa <- readRDS("clean/5_nmr_feces_wide_dm.RDS") %>%
  dplyr::rename(sampleno = sampleid) %>%
  mutate(sampleno = str_extract(sampleno, "[0-9]{2}$")) %>% #prepare for correlation
  dplyr::rename_with(~paste0("m_",.x), !starts_with("sampleno"))

# filter and normalize metagenomic data

g_species_reads_filtered <- breport_reads_filtered %>%
  filter(rank == "S") %>%
  dplyr::select(-rank)

g_il_species_reads <- g_species_reads_filtered %>%
  filter_ileum() %>%
  prepare_omics_for_correlation(abundance_column = "reads") %>%
  dplyr::rename_with(~paste0("g_",.x), !starts_with("sampleno"))

g_fa_species_reads <- g_species_reads_filtered %>%
  filter_faeces() %>%
  prepare_omics_for_correlation(abundance_column = "reads") %>%
  dplyr::rename_with(~paste0("g_",.x), !starts_with("sampleno"))

g_il_kegg_reads <- kegg_reads_filtered %>%
  filter_ileum() %>%
  mutate(reads = as.integer(reads)) %>%
  prepare_omics_for_correlation(abundance_column = "reads") %>%
  dplyr::rename_with(~paste0("g_",.x), !starts_with("sampleno"))

g_fa_kegg_reads <- kegg_reads_filtered %>%
  filter_faeces() %>%
  mutate(reads = as.integer(reads)) %>%
  prepare_omics_for_correlation(abundance_column = "reads") %>%
  dplyr::rename_with(~paste0("g_",.x), !starts_with("sampleno"))

# split datasets in host, micro and pea

proteins_norm_imp_log2_pea <- proteins_norm_imp_log2 %>%
  filter(origin %in% c("pea")) %>%
  separate(proteinid, into = c("p1", "p2", "p3"), sep = "\\|") %>%
  dplyr::select(-p1, -p2, proteinid = p3) # to shorten the final name

# split datasets for ileum and faeces

p_il_pea <- proteins_norm_imp_log2_pea %>%
  filter_ileum() %>%
  prepare_omics_for_correlation(abundance_column = "log2")# %>%
  #dplyr::rename_with(~paste0("pea_",.x), !starts_with("sampleno"))

p_fa_pea <- proteins_norm_imp_log2_pea %>%
  filter_faeces() %>%
  prepare_omics_for_correlation(abundance_column = "log2") #%>%
  #dplyr::rename_with(~paste0("pea_",.x), !starts_with("sampleno"))


p_il_host_kegg <- host_kegg_norm_imp_log2 %>%
  filter_ileum() %>%
  prepare_omics_for_correlation(abundance_column = "log2") %>%
  dplyr::rename_with(~paste0("ssc_",.x), !starts_with("sampleno"))

p_fa_host_kegg <- host_kegg_norm_imp_log2 %>%
  filter_faeces() %>%
  prepare_omics_for_correlation(abundance_column = "log2") %>%
  dplyr::rename_with(~paste0("ssc_",.x), !starts_with("sampleno"))


p_il_kegg <- kegg_norm_imp_log2 %>%
  filter_ileum() %>%
  prepare_omics_for_correlation(abundance_column = "log2") %>%
  dplyr::rename_with(~paste0("p_",.x), !starts_with("sampleno"))

p_fa_kegg <- kegg_norm_imp_log2 %>%
  filter_faeces() %>%
  prepare_omics_for_correlation(abundance_column = "log2") %>%
  dplyr::rename_with(~paste0("p_",.x), !starts_with("sampleno"))

p_genus_log2 <- taxonomy_norm_imp_intensity %>%
  filter(rank == "G") %>% 
  mutate(log2 = log2(intensity+1)) %>%
  dplyr::select(-rank, -intensity)

p_il_genus <- p_genus_log2 %>%
  filter_ileum() %>%
  prepare_omics_for_correlation(abundance_column = "log2") %>%
  dplyr::rename_with(~paste0("p_",.x), !starts_with("sampleno"))

p_fa_genus <- p_genus_log2 %>%
  filter_faeces() %>%
  prepare_omics_for_correlation(abundance_column = "log2") %>%
  dplyr::rename_with(~paste0("p_",.x), !starts_with("sampleno"))

# create matrices

sampleno_intersect_il <- Reduce(intersect, list(pcd$sampleno, 
                                             digesta_analysis$sampleno, 
                                             g_il_species_reads$sampleno, 
                                             g_il_kegg_reads$sampleno,
                                             p_il_genus$sampleno,
                                             p_il_kegg$sampleno,
                                             p_il_host_kegg$sampleno,
                                             p_il_pea$sampleno,
                                             nmr_il$sampleno,
                                             enzymes_il$sampleno))

sampleno_intersect_fa <- Reduce(intersect, list(hindgut$sampleno, 
                                                faeces_analysis$sampleno, 
                                                g_fa_species_reads$sampleno, 
                                                g_fa_kegg_reads$sampleno,
                                                p_fa_genus$sampleno,
                                                p_fa_kegg$sampleno,
                                                p_fa_host_kegg$sampleno,
                                                p_fa_pea$sampleno,
                                                nmr_fa$sampleno,
                                                enzymes_fa$sampleno))

# for ileum
n_pcd_mat <- prepare_for_mixomics(pcd, sampleno_intersect_il)
n_dig_mat <- prepare_for_mixomics(digesta_analysis, sampleno_intersect_il)
g_il_genus_mat <- prepare_for_mixomics(g_il_species_reads, sampleno_intersect_il) %>% 
  t() %>%
  rlog() %>%
  t()
g_il_kegg_mat <- prepare_for_mixomics(g_il_kegg_reads, sampleno_intersect_il) %>% 
  t() %>% 
  rlog() %>%
  t()
p_il_genus_mat <- prepare_for_mixomics(p_il_genus, sampleno_intersect_il)
p_il_kegg_mat <- prepare_for_mixomics(p_il_kegg, sampleno_intersect_il)
p_il_host_mat <- prepare_for_mixomics(p_il_host_kegg, sampleno_intersect_il)
p_il_pea_mat <- prepare_for_mixomics(p_il_pea, sampleno_intersect_il)
m_il_mat <- prepare_for_mixomics(nmr_il, sampleno_intersect_il)
enz_il_mat <- prepare_for_mixomics(enzymes_il, sampleno_intersect_il)

# for faeces
n_hg_mat <- prepare_for_mixomics(hindgut, sampleno_intersect_fa)
n_fa_mat <- prepare_for_mixomics(faeces_analysis, sampleno_intersect_fa)
g_fa_genus_mat <- prepare_for_mixomics(g_fa_species_reads, sampleno_intersect_fa) %>% 
  t() %>%
  rlog() %>%
  t()
g_fa_kegg_mat <- prepare_for_mixomics(g_fa_kegg_reads, sampleno_intersect_fa) %>% 
  t() %>% 
  rlog() %>%
  t()
p_fa_genus_mat <- prepare_for_mixomics(p_fa_genus, sampleno_intersect_fa)
p_fa_kegg_mat <- prepare_for_mixomics(p_fa_kegg, sampleno_intersect_fa)
p_fa_host_mat <- prepare_for_mixomics(p_fa_host_kegg, sampleno_intersect_fa)
p_fa_pea_mat <- prepare_for_mixomics(p_fa_pea, sampleno_intersect_fa)
m_fa_mat <- prepare_for_mixomics(nmr_fa, sampleno_intersect_fa)
enz_fa_mat <- prepare_for_mixomics(enzymes_fa, sampleno_intersect_fa)

# create final matrices

cov_il <- meta_il %>%
  dplyr::select(animal, period)

n_il_micro <- dataPreprocess(cbind(n_pcd_mat, n_dig_mat), covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
n_il_host <- dataPreprocess(cbind(n_pcd_mat, n_dig_mat, enz_il_mat), covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
g_il <- dataPreprocess(cbind(g_il_genus_mat, g_il_kegg_mat), covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
p_il_micro <- dataPreprocess(cbind(p_il_genus_mat, p_il_kegg_mat), covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
p_il_host <- dataPreprocess(p_il_host_mat, covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
p_il_pea <- dataPreprocess(p_il_pea_mat, covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
m_il <- dataPreprocess(m_il_mat, covariates = cov_il, is_cv = TRUE, cv_quantile = .1)
#enz_il <- dataPreprocess(enz_il_mat, covariates = cov_il, is_cv = TRUE, cv_quantile = .0)

cov_fa <- meta_fa %>%
  dplyr::select(animal, period)

n_fa_micro <- dataPreprocess(cbind(n_hg_mat, n_fa_mat), covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
n_fa_host <- dataPreprocess(cbind(n_hg_mat, n_fa_mat, enz_fa_mat), covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
g_fa <- dataPreprocess(cbind(g_fa_genus_mat, g_fa_kegg_mat), covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
p_fa_micro <- dataPreprocess(cbind(p_fa_genus_mat, p_fa_kegg_mat), covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
p_fa_host <- dataPreprocess(p_fa_host_mat, covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
p_fa_pea <- dataPreprocess(p_fa_pea_mat, covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
m_fa <- dataPreprocess(m_fa_mat, covariates = cov_fa, is_cv = TRUE, cv_quantile = .1)
#enz_fa <- dataPreprocess(enz_fa_mat, covariates = cov_fa, is_cv = TRUE, cv_quantile = .0)

# create list objects



il_combined <- list("Nutrition" = n_il_host,
                    "Metagenomics" = g_il,
                    "Metaproteomics" = p_il_micro,
                    "Host" = p_il_host,
                    "Pea" = p_il_pea,
                    "Metabolomics" = m_il)

saveRDS(il_combined, "temp/70_il_combined_list.RDS")



fa_combined <- list("Nutrition" = n_fa_host,
                    "Metagenomics" = g_fa,
                    "Metaproteomics" = p_fa_micro,
                    "Host" = p_fa_host,
                    "Pea" = p_fa_pea,
                    "Metabolomics" = m_fa)

saveRDS(fa_combined, "temp/70_fa_combined_list.RDS")