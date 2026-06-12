
library(here)

source(here("20_nutrition_functions.R"))
source("E:/R/source/ggplot2_theme.R")

# load data

performance <- read_in_nutrition(here("data/nutrition_performance.txt")) 

pcd <- read_in_nutrition(here("data/nutrition_pcd.txt"))

pcd_g <- read_in_nutrition(here("data/nutrition_pcd_g.txt"))

ttd <- read_in_nutrition(here("data/nutrition_ttd.txt")) 

ttd_g <- read_in_nutrition(here("data/nutrition_ttd_g.txt"))

hindgut <- read_in_nutrition(here("data/nutrition_hindgut.txt")) 

hindgut_g <- read_in_nutrition(here("data/nutrition_hindgut_g.txt"))

retention <- read_in_nutrition(here("data/nutrition_retention.txt"))

excretion <- read_in_nutrition(here("data/nutrition_excreted.txt"))

# saveswitch

save = FALSE
# save = TRUE

# performance data
# LMZ

performance_lmz <- combined_comparison_and_results(performance, selected_response = "lmz", 
                              transformation = "log", type = "s", response_name = "weight gain",
                              y_axis = "g/d", save_name = "performance_lmz", save = save)

# precaecal digestibility

pcd_dm <- combined_comparison_and_results(pcd, selected_response = "dm", 
                                                   transformation = "none", type = "s", response_name = "DM",
                                                   y_axis = "digestibility [%]", save_name = "pcd_dm", save = save)

pcd_starch <- combined_comparison_and_results(pcd, selected_response = "total_starch", 
                                          transformation = "tukey", type = "s", response_name = "total starch",
                                          y_axis = "digestibility [%]", save_name = "pcd_starch", save = save)

pcd_tdf <- combined_comparison_and_results(pcd, selected_response = "tdf", 
                                              transformation = "none", type = "s", response_name = "TDF",
                                              y_axis = "digestibility [%]", save_name = "pcd_tdf", save = save)

pcd_cp <- combined_comparison_and_results(pcd, selected_response = "cp", 
                                              transformation = "tukey", type = "hh", response_name = "CP",
                                              y_axis = "digestibility [%]", save_name = "pcd_cp", save = save)

pcd_n <- combined_comparison_and_results(pcd, selected_response = "n", 
                                          transformation = "tukey", type = "hh", response_name = "N",
                                          y_axis = "digestibility [%]", save_name = "pcd_n", save = save)

pcd_insp6 <- combined_comparison_and_results(pcd, selected_response = "insp6", 
                                         transformation = "none", type = "s", response_name = "InsP6",
                                         y_axis = "digestibility [%]", save_name = "pcd_insp6", save = save)

pcd_p <- combined_comparison_and_results(pcd, selected_response = "p", 
                                          transformation = "tukey", type = "s", response_name = "P",
                                          y_axis = "digestibility [%]", save_name = "pcd_p", save = save)

pcd_ca <- combined_comparison_and_results(pcd, selected_response = "ca", 
                                         transformation = "tukey", type = "s", response_name = "Ca",
                                         y_axis = "digestibility [%]", save_name = "pcd_ca", save = save)

pcd_ge <- combined_comparison_and_results(pcd, selected_response = "ge", 
                                          transformation = "none", type = "s", response_name = "GE",
                                          y_axis = "digestibility [%]", save_name = "pcd_ge", save = save)

pcd_arg <- combined_comparison_and_results(pcd, selected_response = "arg", 
                                          transformation = "tukey", type = "s", response_name = "Arg",
                                          y_axis = "digestibility [%]", save_name = "pcd_arg", save = save)

pcd_his <- combined_comparison_and_results(pcd, selected_response = "his", 
                                           transformation = "tukey", type = "s", response_name = "His",
                                           y_axis = "digestibility [%]", save_name = "pcd_his", save = save)

pcd_ile <- combined_comparison_and_results(pcd, selected_response = "ile", 
                                           transformation = "none", type = "s", response_name = "Ile",
                                           y_axis = "digestibility [%]", save_name = "pcd_ile", save = save)

pcd_leu <- combined_comparison_and_results(pcd, selected_response = "leu", 
                                           transformation = "tukey", type = "s", response_name = "Leu",
                                           y_axis = "digestibility [%]", save_name = "pcd_leu", save = save)

pcd_lys <- combined_comparison_and_results(pcd, selected_response = "lys", 
                                           transformation = "tukey", type = "s", response_name = "Lys",
                                           y_axis = "digestibility [%]", save_name = "pcd_lys", save = save)

pcd_met <- combined_comparison_and_results(pcd, selected_response = "met", 
                                           transformation = "none", type = "s", response_name = "Met",
                                           y_axis = "digestibility [%]", save_name = "pcd_met", save = save)

pcd_phe <- combined_comparison_and_results(pcd, selected_response = "phe", 
                                           transformation = "none", type = "s", response_name = "Phe",
                                           y_axis = "digestibility [%]", save_name = "pcd_phe", save = save)

pcd_thr <- combined_comparison_and_results(pcd, selected_response = "thr", 
                                           transformation = "none", type = "s", response_name = "Thr",
                                           y_axis = "digestibility [%]", save_name = "pcd_thr", save = save)

pcd_val <- combined_comparison_and_results(pcd, selected_response = "val", 
                                           transformation = "none", type = "s", response_name = "Val",
                                           y_axis = "digestibility [%]", save_name = "pcd_val", save = save)

pcd_ala <- combined_comparison_and_results(pcd, selected_response = "ala", 
                                           transformation = "tukey2", type = "s", response_name = "Ala",
                                           y_axis = "digestibility [%]", save_name = "pcd_ala", save = save)

pcd_asp <- combined_comparison_and_results(pcd, selected_response = "asp", 
                                           transformation = "tukey2", type = "s", response_name = "Asp",
                                           y_axis = "digestibility [%]", save_name = "pcd_asp", save = save)

pcd_cys <- combined_comparison_and_results(pcd, selected_response = "cys", 
                                           transformation = "tukey", type = "s", response_name = "Cys",
                                           y_axis = "digestibility [%]", save_name = "pcd_cys", save = save)

pcd_glu <- combined_comparison_and_results(pcd, selected_response = "glu", 
                                           transformation = "none", type = "s", response_name = "Glu",
                                           y_axis = "digestibility [%]", save_name = "pcd_asp", save = save)

pcd_gly <- combined_comparison_and_results(pcd, selected_response = "gly", 
                                           transformation = "none", type = "s", response_name = "Gly",
                                           y_axis = "digestibility [%]", save_name = "pcd_gly", save = save)

pcd_pro <- combined_comparison_and_results(pcd, selected_response = "pro", 
                                           transformation = "gauss", type = "s", response_name = "Pro",
                                           y_axis = "digestibility [%]", save_name = "pcd_pro", save = save)

pcd_ser <- combined_comparison_and_results(pcd, selected_response = "ser", 
                                           transformation = "tukey", type = "s", response_name = "Ser",
                                           y_axis = "digestibility [%]", save_name = "pcd_ser", save = save)

pcd_tyr <- combined_comparison_and_results(pcd, selected_response = "tyr", 
                                           transformation = "tukey", type = "s", response_name = "Tyr",
                                           y_axis = "digestibility [%]", save_name = "pcd_tyr", save = save)

# create AA digestibility plot

table_aa <- cbind(pcd_arg, pcd_his[2], pcd_ile[2], pcd_leu[2], pcd_lys[2], pcd_met[2], pcd_phe[2], pcd_thr[2], pcd_val[2], pcd_ala[2], pcd_asp[2], pcd_cys[2], pcd_glu[2], pcd_gly[2], pcd_pro[2], pcd_ser[2], pcd_tyr[2]) %>%
  pivot_longer(-diet, names_to = "aa", values_to = "digestibility") %>%
  filter(diet != "P-value") %>%
  pivot_wider(names_from="diet", values_from="digestibility") %>%
  dplyr::rename("sem" = `Pooled SEM`) %>%
  pivot_longer(-c("aa", "sem"), names_to = "diet", values_to = "digestibility") %>%
  mutate(letters = str_extract(digestibility, "[a-z]+"),
         digestibility = as.numeric(str_remove(digestibility, "[a-z]+")),
         sem = as.numeric(sem))


ggplot(table_aa, aes(x = factor(aa, levels = c("Arg", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Thr", "Val", "Ala", "Asp", "Cys", "Glu", "Gly", "Pro", "Ser", "Tyr")), y = digestibility, fill = diet)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = digestibility + sem, ymin = digestibility - sem),
                width = 0, size = 0.7, position = position_dodge(width = 0.9)) +
  geom_text(aes(y = (digestibility + sem + 2), label = str_trim(letters)), size = 2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("black", "grey20", "grey40", "grey60")) +
  scale_y_continuous(limits = c(0, 1.1*max(table_aa$digestibility)), expand = expansion(mult = c(0, .1))) +
  labs(x = "", y = "Digestibility", title = "pc Amino acid digestibility", fill = "Diet")

save_big("21_pcd_aa_digestibility")
  

# precaecal digestibility g/kg

pcd_dm_g <- combined_comparison_and_results(pcd_g, selected_response = "dm", 
                                          transformation = "test", type = "s", response_name = "DM",
                                          y_axis = "digested [g/kg DM]", save_name = "pcd_dm_g", save = save)

pcd_starch_g <- combined_comparison_and_results(pcd_g, selected_response = "total_starch", 
                                              transformation = "test", type = "s", response_name = "total starch",
                                              y_axis = "digested [g/kg DM]", save_name = "pcd_starch_g", save = save)

pcd_tdf_g <- combined_comparison_and_results(pcd_g, selected_response = "tdf", 
                                                transformation = "none", type = "s", response_name = "TDF",
                                                y_axis = "digested [g/kg DM]", save_name = "pcd_tdf_g", save = save)

pcd_cp_g <- combined_comparison_and_results(pcd_g, selected_response = "cp", 
                                          transformation = "test", type = "hh", response_name = "CP",
                                          y_axis = "digested [g/kg DM]", save_name = "pcd_cp_g", save = save)

pcd_n_g <- combined_comparison_and_results(pcd_g, selected_response = "n", 
                                         transformation = "test", type = "hh", response_name = "N",
                                         y_axis = "digested [g/kg DM]", save_name = "pcd_n_g", save = save)

pcd_insp6p_g <- combined_comparison_and_results(pcd_g, selected_response = "insp6_p", 
                                           transformation = "none", type = "hh", response_name = "InsP6-P",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_insp6p_g", save = save)

pcd_p_g <- combined_comparison_and_results(pcd_g, selected_response = "p", 
                                         transformation = "test", type = "s", response_name = "P",
                                         y_axis = "digested [g/kg DM]", save_name = "pcd_p_g", save = save)

pcd_ca_g <- combined_comparison_and_results(pcd_g, selected_response = "ca", 
                                          transformation = "test", type = "s", response_name = "Ca",
                                          y_axis = "digested [g/kg DM]", save_name = "pcd_ca_g", save = save)

pcd_ge_g <- combined_comparison_and_results(pcd_g, selected_response = "ge", 
                                          transformation = "test", type = "s", response_name = "GE",
                                          y_axis = "digested [MJ/kg DM]", save_name = "pcd_ge_g", save = save)

pcd_arg_g <- combined_comparison_and_results(pcd_g, selected_response = "arg", 
                                           transformation = "test", type = "s", response_name = "Arg",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_arg_g", save = save)

pcd_his_g <- combined_comparison_and_results(pcd_g, selected_response = "his", 
                                           transformation = "test", type = "s", response_name = "His",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_his_g", save = save)

pcd_ile_g <- combined_comparison_and_results(pcd_g, selected_response = "ile", 
                                           transformation = "test", type = "s", response_name = "Ile",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_ile_g", save = save)

pcd_leu_g <- combined_comparison_and_results(pcd_g, selected_response = "leu", 
                                           transformation = "test", type = "s", response_name = "Leu",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_leu_g", save = save)

pcd_lys_g <- combined_comparison_and_results(pcd_g, selected_response = "lys", 
                                           transformation = "test", type = "s", response_name = "Lys",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_lys_g", save = save)

pcd_met_g <- combined_comparison_and_results(pcd_g, selected_response = "met", 
                                           transformation = "test", type = "s", response_name = "Met",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_met_g", save = save)

pcd_phe_g <- combined_comparison_and_results(pcd_g, selected_response = "phe", 
                                           transformation = "test", type = "s", response_name = "Phe",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_phe_g", save = save)

pcd_thr_g <- combined_comparison_and_results(pcd_g, selected_response = "thr", 
                                           transformation = "test", type = "s", response_name = "Thr",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_thr_g", save = save)

pcd_val_g <- combined_comparison_and_results(pcd_g, selected_response = "val", 
                                           transformation = "test", type = "s", response_name = "Val",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_val_g", save = save)

pcd_ala_g <- combined_comparison_and_results(pcd_g, selected_response = "ala", 
                                           transformation = "test", type = "s", response_name = "Ala",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_ala_g", save = save)

pcd_asp_g <- combined_comparison_and_results(pcd_g, selected_response = "asp", 
                                           transformation = "test", type = "s", response_name = "Asp",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_asp_g", save = save)

pcd_cys_g <- combined_comparison_and_results(pcd_g, selected_response = "cys", 
                                           transformation = "test", type = "s", response_name = "Cys",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_cys_g", save = save)

pcd_glu_g <- combined_comparison_and_results(pcd_g, selected_response = "glu", 
                                           transformation = "test", type = "s", response_name = "Glu",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_asp_g", save = save)

pcd_gly_g <- combined_comparison_and_results(pcd_g, selected_response = "gly", 
                                           transformation = "test", type = "s", response_name = "Gly",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_gly_g", save = save)

pcd_pro_g <- combined_comparison_and_results(pcd_g, selected_response = "pro", 
                                           transformation = "test", type = "s", response_name = "Pro",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_pro_g", save = save)

pcd_ser_g <- combined_comparison_and_results(pcd_g, selected_response = "ser", 
                                           transformation = "test", type = "s", response_name = "Ser",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_ser_g", save = save)

pcd_tyr_g <- combined_comparison_and_results(pcd_g, selected_response = "tyr", 
                                           transformation = "test", type = "s", response_name = "Tyr",
                                           y_axis = "digested [g/kg DM]", save_name = "pcd_tyr_g", save = save)

# total tract digestibility

ttd_dm <- combined_comparison_and_results(ttd, selected_response = "dm", 
                                           transformation = "none", type = "s", response_name = "DM",
                                           y_axis = "digestibility [%]", save_name = "ttd_dm", save = save)

ttd_tdf <- combined_comparison_and_results(ttd, selected_response = "tdf", 
                                          transformation = "none", type = "s", response_name = "TDF",
                                          y_axis = "digestibility [%]", save_name = "ttd_tdf", save = save)

ttd_cp <- combined_comparison_and_results(ttd, selected_response = "cp", 
                                          transformation = "tukey2", type = "hh", response_name = "CP",
                                          y_axis = "digestibility [%]", save_name = "ttd_cp", save = save)

ttd_n <- combined_comparison_and_results(ttd, selected_response = "n", 
                                          transformation = "tukey2", type = "hh", response_name = "N",
                                          y_axis = "digestibility [%]", save_name = "ttd_n", save = save)

ttd_insp6 <- combined_comparison_and_results(ttd, selected_response = "insp6", 
                                         transformation = "gauss", type = "s", response_name = "InsP6",
                                         y_axis = "digestibility [%]", save_name = "ttd_insp6", save = save)

ttd_p <- combined_comparison_and_results(ttd, selected_response = "p", 
                                          transformation = "none", type = "s", response_name = "P",
                                          y_axis = "digestibility [%]", save_name = "ttd_p", save = save)

ttd_ca <- combined_comparison_and_results(ttd, selected_response = "ca", 
                                         transformation = "tukey", type = "s", response_name = "Ca",
                                         y_axis = "digestibility [%]", save_name = "ttd_ca", save = save)

ttd_k <- combined_comparison_and_results(ttd, selected_response = "k", 
                                          transformation = "none", type = "s", response_name = "K",
                                          y_axis = "digestibility [%]", save_name = "ttd_k", save = save)

ttd_ge <- combined_comparison_and_results(ttd, selected_response = "ge", 
                                         transformation = "none", type = "s", response_name = "GE",
                                         y_axis = "digestibility [%]", save_name = "ttd_ge", save = save)

# ttd digestibility g/kg

ttd_dm_g <- combined_comparison_and_results(ttd_g, selected_response = "dm", 
                                          transformation = "test", type = "s", response_name = "DM",
                                          y_axis = "digested [g/kg DM]", save_name = "ttd_dm_g", save = save)

ttd_tdf_g <- combined_comparison_and_results(ttd_g, selected_response = "tdf", 
                                            transformation = "none", type = "s", response_name = "TDF",
                                            y_axis = "digested [g/kg DM]", save_name = "ttd_tdf_g", save = save)

ttd_cp_g <- combined_comparison_and_results(ttd_g, selected_response = "cp", 
                                          transformation = "test", type = "hh", response_name = "CP",
                                          y_axis = "digested [g/kg DM]", save_name = "ttd_cp_g", save = save)

ttd_n_g <- combined_comparison_and_results(ttd_g, selected_response = "n", 
                                         transformation = "test", type = "hh", response_name = "N",
                                         y_axis = "digested [g/kg DM]", save_name = "ttd_n_g", save = save)

ttd_insp6p_g <- combined_comparison_and_results(ttd_g, selected_response = "insp6_p", 
                                           transformation = "gauss", type = "s", response_name = "InsP6-P",
                                           y_axis = "digested [g/kg DM]", save_name = "ttd_insp6p_g", save = save)

ttd_p_g <- combined_comparison_and_results(ttd_g, selected_response = "p", 
                                         transformation = "test", type = "s", response_name = "P",
                                         y_axis = "digested [g/kg DM]", save_name = "ttd_p_g", save = save)

ttd_ca_g <- combined_comparison_and_results(ttd_g, selected_response = "ca", 
                                          transformation = "test", type = "s", response_name = "Ca",
                                          y_axis = "digested [g/kg DM]", save_name = "ttd_ca_g", save = save)

ttd_k_g <- combined_comparison_and_results(ttd_g, selected_response = "k", 
                                         transformation = "test", type = "s", response_name = "K",
                                         y_axis = "digested [g/kg DM]", save_name = "ttd_k_g", save = save)

ttd_ge_g <- combined_comparison_and_results(ttd_g, selected_response = "ge", 
                                          transformation = "test", type = "s", response_name = "GE",
                                          y_axis = "digested [MJ/kg DM]", save_name = "ttd_ge_g", save = save)

# hindgut

hindgut_dm <- combined_comparison_and_results(hindgut, selected_response = "dm", 
                                              transformation = "test", type = "s", response_name = "DM",
                                              y_axis = "disappearance [%]", save_name = "hindgut_dm", save = save)

hindgut_tdf <- combined_comparison_and_results(hindgut, selected_response = "tdf", 
                                              transformation = "tukey2", type = "s", response_name = "TDF",
                                              y_axis = "disappearance [%]", save_name = "hindgut_tdf", save = save)

hindgut_cp <- combined_comparison_and_results(hindgut, selected_response = "cp", 
                                              transformation = "none", type = "s", response_name = "CP",
                                              y_axis = "disappearance [%]", save_name = "hindgut_cp", save = save)

hindgut_n <- combined_comparison_and_results(hindgut, selected_response = "n", 
                                              transformation = "none", type = "s", response_name = "N",
                                              y_axis = "disappearance [%]", save_name = "hindgut_n", save = save)

hindgut_insp6 <- combined_comparison_and_results(hindgut, selected_response = "insp6", 
                                             transformation = "gauss", type = "s", response_name = "InsP6",
                                             y_axis = "disappearance [%]", save_name = "hindgut_insp6", save = save)

hindgut_p <- combined_comparison_and_results(hindgut, selected_response = "p", 
                                              transformation = "gauss", type = "s", response_name = "P",
                                              y_axis = "disappearance [%]", save_name = "hindgut_p", save = save)

hindgut_ca <- combined_comparison_and_results(hindgut, selected_response = "ca", 
                                             transformation = "none", type = "s", response_name = "Ca",
                                             y_axis = "disappearance [%]", save_name = "hindgut_ca", save = save)

hindgut_ge <- combined_comparison_and_results(hindgut, selected_response = "ge", 
                                              transformation = "none", type = "s", response_name = "GE",
                                              y_axis = "disappearance [%]", save_name = "hindgut_ge", save = save)

# hindgut g/kg DM

hindgut_dm_g <- combined_comparison_and_results(hindgut_g, selected_response = "dm", 
                                              transformation = "tukey2", type = "s", response_name = "DM",
                                              y_axis = "disappearance [g/kg DM]", save_name = "hindgut_dm_g", save = save)

hindgut_tdf_g <- combined_comparison_and_results(hindgut_g, selected_response = "tdf", 
                                                transformation = "tukey2", type = "s", response_name = "TDF",
                                                y_axis = "disappearance [g/kg DM]", save_name = "hindgut_tdf_g", save = save)

hindgut_cp_g <- combined_comparison_and_results(hindgut_g, selected_response = "cp", 
                                              transformation = "tukey2", type = "s", response_name = "CP",
                                              y_axis = "disappearance [g/kg DM]", save_name = "hindgut_cp_g", save = save)

hindgut_n_g <- combined_comparison_and_results(hindgut_g, selected_response = "n", 
                                             transformation = "tukey2", type = "s", response_name = "N",
                                             y_axis = "disappearance [g/kg DM]", save_name = "hindgut_n_g", save = save)

hindgut_insp6p_g <- combined_comparison_and_results(hindgut_g, selected_response = "insp6_p", 
                                               transformation = "gauss", type = "s", response_name = "InsP6-P",
                                               y_axis = "disappearance [g/kg DM]", save_name = "hindgut_insp6p_g", save = save)

hindgut_p_g <- combined_comparison_and_results(hindgut_g, selected_response = "p", 
                                             transformation = "gauss", type = "s", response_name = "P",
                                             y_axis = "disappearance [g/kg DM]", save_name = "hindgut_p_g", save = save)

hindgut_ca_g <- combined_comparison_and_results(hindgut_g, selected_response = "ca", 
                                              transformation = "none", type = "s", response_name = "Ca",
                                              y_axis = "disappearance [g/kg DM]", save_name = "hindgut_ca_g", save = save)

hindgut_ge_g <- combined_comparison_and_results(hindgut_g, selected_response = "ge", 
                                              transformation = "tukey2", type = "s", response_name = "GE",
                                              y_axis = "disappearance [MJ/kg DM]", save_name = "hindgut_ge_g", save = save)

# retention

retention_p_g_d <- combined_comparison_and_results(retention, selected_response = "p_g_d", 
                                         transformation = "tukey2", type = "s", response_name = "P g/d",
                                         y_axis = "retention [g/d]", save_name = "retention_p_g_d", save = save)

retention_p <- combined_comparison_and_results(retention, selected_response = "p_%", 
                                         transformation = "tukey", type = "s", response_name = "P",
                                         y_axis = "retention [%]", save_name = "retention_p", save = save)

retention_p_g_kg <- combined_comparison_and_results(retention, selected_response = "p_g_kg", 
                                                   transformation = "tukey2", type = "s", response_name = "P g/kg",
                                                   y_axis = "retention [g/kg DM]", save_name = "retention_p_g_kg", save = save)

retention_ca_g_d <- combined_comparison_and_results(retention, selected_response = "ca_g_d", 
                                     transformation = "tukey2", type = "s", response_name = "Ca g/d",
                                     y_axis = "retention [g/d]", save_name = "retention_ca_g_d", save = save)

retention_ca <- combined_comparison_and_results(retention, selected_response = "ca_%", 
                                                transformation = "tukey2", type = "s", response_name = "Ca",
                                                y_axis = "retention [%]", save_name = "retention_ca", save = save)

retention_ca_g_kg <- combined_comparison_and_results(retention, selected_response = "ca_g_kg", 
                                                transformation = "tukey2", type = "s", response_name = "Ca g/kg",
                                                y_axis = "retention [g/kg DM]", save_name = "retention_ca_g_kg", save = save)

retention_n_g_d <- combined_comparison_and_results(retention[-32,], selected_response = "n_g_d", 
                                                     transformation = "test", type = "s", response_name = "N g/d",
                                                     y_axis = "retention [g/d]", save_name = "retention_n_g_d", save = save)

retention_n_percent <- combined_comparison_and_results(retention[-32,], selected_response = "n_%", 
                                                   transformation = "test", type = "s", response_name = "N %",
                                                   y_axis = "retention [%]", save_name = "retention_n_percent", save = save)

# excreted minerals

excretion_p_mg_d <- combined_comparison_and_results(excretion, selected_response = "p_mg_d", 
                                                       transformation = "test", type = "none", response_name = "P mg/d",
                                                       y_axis = "excretion [mg/d]", save_name = "excretion_p_mg_d", save = save)

excretion_p_mg_kg_bw_d <- combined_comparison_and_results(excretion, selected_response = "p_mg_kg_bw_d", 
                                                    transformation = "test", type = "none", response_name = "P mg/kg BW/d",
                                                    y_axis = "excretion [mg/kg BW/d]", save_name = "excretion_p_mg_kg_bw_d", save = save)

excretion_ca_mg_d <- combined_comparison_and_results(excretion, selected_response = "ca_mg_d", 
                                                    transformation = "test", type = "none", response_name = "Ca mg/d",
                                                    y_axis = "excretion [mg/d]", save_name = "excretion_ca_mg_d", save = save)

excretion_ca_mg_kg_bw_d <- combined_comparison_and_results(excretion, selected_response = "ca_mg_kg_bw_d", 
                                                          transformation = "test", type = "none", response_name = "Ca mg/kg BW/d",
                                                          y_axis = "excretion [mg/kg BW/d]", save_name = "excretion_ca_mg_kg_bw_d", save = save)

excretion_n_g_d <- combined_comparison_and_results(excretion, selected_response = "n_g_d", 
                                                     transformation = "test", type = "none", response_name = "N g/d",
                                                     y_axis = "excretion [g/d]", save_name = "excretion_n_g_d", save = save)

excretion_urea_g_d <- combined_comparison_and_results(excretion, selected_response = "urea_g_d", 
                                                   transformation = "test", type = "none", response_name = "Urea g/d",
                                                   y_axis = "excretion [g/d]", save_name = "excretion_urea_g_d", save = save)

excretion_urea_n_g_d <- combined_comparison_and_results(excretion, selected_response = "urea_n_g_d", 
                                                      transformation = "test", type = "none", response_name = "Urea-N g/d",
                                                      y_axis = "excretion [g/d]", save_name = "excretion_urea_n_g_d", save = save)

excretion_urea_of_total_n <- combined_comparison_and_results(excretion, selected_response = "urea_of_total_n", 
                                                      transformation = "none", type = "none", response_name = "Urea % of total N",
                                                      y_axis = "excretion [%]", save_name = "excretion_urea_of_total_n", save = save)

excretion_myo_inositol_mg_d <- combined_comparison_and_results(excretion[-32,], selected_response = "myo_inositol_mg_d", 
                                                     transformation = "test", type = "none", response_name = "myo-Inositol mg/d",
                                                     y_axis = "excretion [mg/d]", save_name = "excretion_myo_inositol_mg_d", save = save)

excretion_myo_inositol_mg_kg_bw_d <- combined_comparison_and_results(excretion[-32,], selected_response = "myo_inositol_mg_kg_bw_d", 
                                                               transformation = "test", type = "none", response_name = "myo-Inositol mg/kg BW/d",
                                                               y_axis = "excretion [mg/kg BW/d]", save_name = "excretion_myo_inositol_mg_kg_bw_d", save = save)

# combine tables

performance_table <- cbind(performance_lmz)
write_tsv(performance_table, "tables/21_performance.txt")

pcd_table <- cbind(pcd_dm, pcd_starch[2], pcd_tdf[2], pcd_cp[2], pcd_n[2], pcd_insp6[2], pcd_p[2], pcd_ca[2], pcd_ge[2], pcd_arg[2], pcd_his[2], pcd_ile[2], pcd_leu[2], pcd_lys[2], pcd_met[2], pcd_phe[2], pcd_thr[2], pcd_val[2], pcd_ala[2], pcd_asp[2], pcd_cys[2], pcd_glu[2], pcd_gly[2], pcd_pro[2], pcd_ser[2], pcd_tyr[2])
write_tsv(pcd_table, "tables/21_pcd.txt")

pcd_g_table <- cbind(pcd_dm_g, pcd_starch_g[2], pcd_tdf_g[2], pcd_cp_g[2], pcd_n_g[2], pcd_insp6p_g[2], pcd_p_g[2], pcd_ca_g[2], pcd_ge_g[2], pcd_arg_g[2], pcd_his_g[2], pcd_ile_g[2], pcd_leu_g[2], pcd_lys_g[2], pcd_met_g[2], pcd_phe_g[2], pcd_thr_g[2], pcd_val_g[2], pcd_ala_g[2], pcd_asp_g[2], pcd_cys_g[2], pcd_glu_g[2], pcd_gly_g[2], pcd_pro_g[2], pcd_ser_g[2], pcd_tyr_g[2])
write_tsv(pcd_g_table, "tables/21_pcd_g.txt")

ttd_table <- cbind(ttd_dm, ttd_tdf[2], ttd_cp[2], ttd_n[2], ttd_insp6[2], ttd_p[2], ttd_ca[2], ttd_k[2], ttd_ge[2])
write_tsv(ttd_table, "tables/21_ttd.txt")

ttd_g_table <- cbind(ttd_dm_g, ttd_tdf_g[2], ttd_cp_g[2], ttd_n_g[2], ttd_insp6p_g[2], ttd_p_g[2], ttd_ca_g[2], ttd_k_g[2], ttd_ge_g[2])
write_tsv(ttd_g_table, "tables/21_ttd_g.txt")

hindgut_table <- cbind(hindgut_dm, hindgut_tdf[2], hindgut_cp[2], hindgut_n[2], hindgut_insp6[2], hindgut_p[2], hindgut_ca[2], hindgut_ge[2])
write_tsv(hindgut_table, "tables/21_hindgut.txt")

hindgut_g_table <- cbind(hindgut_dm_g, hindgut_tdf_g[2], hindgut_cp_g[2], hindgut_n_g[2], hindgut_insp6p_g[2], hindgut_p_g[2], hindgut_ca_g[2], hindgut_ge_g[2])
write_tsv(hindgut_g_table, "tables/21_hindgut_g.txt")

retention_table <- cbind(retention_p_g_d, retention_p[2], retention_p_g_kg[2], 
                         retention_ca_g_d[2], retention_ca[2], retention_ca_g_kg[2],
                         retention_n_g_d[2], retention_n_percent[2])
write_tsv(retention_table, "tables/21_retention.txt")

excretion_table <- cbind(excretion_p_mg_d, excretion_p_mg_kg_bw_d[2], excretion_ca_mg_d[2], excretion_ca_mg_kg_bw_d[2], excretion_n_g_d[2], excretion_urea_g_d[2], excretion_urea_n_g_d[2], excretion_urea_of_total_n[2], excretion_myo_inositol_mg_d[2], excretion_myo_inositol_mg_kg_bw_d[2])
write_tsv(excretion_table, "tables/21_excretion.txt")

save_results_table(performance_table, title = "Performance", name = "21_table_performance")
pcd_gt <- save_results_table(pcd_table, title = "pcd digestibility", name = "21_table_pcd")
gtsave(pcd_gt, filename = "21_table_pcd.png", "plots", vwidth = 2000)
pcd_g_gt <- save_results_table(pcd_g_table, title = "pcd digested", name = "21_table_pcd_g")
gtsave(pcd_g_gt, filename = "21_table_pcd_d.png", "plots", vwidth = 2000)
save_results_table(ttd_table, title = "TTD", name = "21_table_ttd")
save_results_table(ttd_g_table, title = "TTD", name = "21_table_ttd_g")
save_results_table(hindgut_table, title = "hindgut disappearance", name = "21_table_hindgut")
save_results_table(hindgut_g_table, title = "hindgut disappearance", name = "21_table_hindgut_g")
save_results_table(retention_table, title = "retention", name = "21_table_retention")
save_results_table(excretion_table, title = "excretion", name = "21_table_excretion")


#save_output()