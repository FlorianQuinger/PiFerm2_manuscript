
set.seed(1112)
library(here)
library(tidyverse)

# import other scripts

source("E:/R/source/ggplot2_theme.R")

# load data

cara_il <- readRDS("clean/12_cara_il.RDS")
cara_fa <- readRDS("clean/12_cara_fa.RDS")
carb_il <- readRDS("clean/11_carb_il.RDS")
carb_fa <- readRDS("clean/11_carb_fa.RDS")
try_il <- readRDS("clean/13_try_il.RDS")
try_fa <- readRDS("clean/13_try_fa.RDS")
chy_il <- readRDS("clean/14_chy_il.RDS")
chy_fa <- readRDS("clean/14_chy_fa.RDS")
amy_il <- readRDS("clean/15_amy_il.RDS")
amy_fa <- readRDS("clean/15_amy_fa.RDS")

# load dm data 

dm <- read_tsv("data/dry_matter.txt") %>%
  mutate(sampleid = as.character(sampleid))

# combine files

cara <- rbind(cara_fa, cara_il) %>%
  dplyr::select(sampleid, cara_fm = mean_iu_g)
carb <- rbind(carb_fa, carb_il) %>%
  dplyr::select(sampleid, carb_fm = mean_iu_g)
try <- rbind(try_fa, try_il) %>%
  dplyr::select(sampleid, try_fm = mean_iu_g)
chy <- rbind(chy_fa, chy_il) %>%
  dplyr::select(sampleid, chy_fm = mean_iu_g)
amy <- rbind(amy_fa, amy_il) %>%
  dplyr::select(sampleid, amy_fm = mean_iu_g)

# combine further

enzymes <- dm %>%
  arrange(sampleid) %>%
  inner_join(cara, by = "sampleid") %>%
  inner_join(carb, by = "sampleid") %>%
  inner_join(try, by = "sampleid") %>%
  inner_join(chy, by = "sampleid") %>%
  inner_join(amy, by = "sampleid")

enzymes_dm <- enzymes %>%
  mutate(cara_dm = cara_fm * dry_matter/100) %>%
  mutate(carb_dm = carb_fm * dry_matter/100) %>%
  mutate(try_dm = try_fm * dry_matter/100) %>%
  mutate(chy_dm = chy_fm * dry_matter/100) %>%
  mutate(amy_dm = amy_fm * dry_matter/100)

saveRDS(enzymes_dm, "clean/16_enzymes_dm.RDS")

# calculate ileal digesta/faeces ratio

enzymes_fa <- enzymes_dm %>%
  filter(str_detect(sampleid, "^3")) %>%
  mutate(sampleid = str_replace(sampleid, "^3", "1"))

enzymes_il <- enzymes_dm %>%
  filter(str_detect(sampleid, "^1"))

enzymes_ratio <- enzymes_il %>%
  inner_join(enzymes_fa, by = "sampleid") %>%
  mutate(cara_ratio = cara_dm.x / cara_dm.y,
         carb_ratio = carb_dm.x / carb_dm.y,
         try_ratio = try_dm.x / try_dm.y,
         chy_ratio = chy_dm.x / chy_dm.y,
         amy_ratio = amy_dm.x / amy_dm.y) %>%
  dplyr::select(sampleid, ends_with("ratio"))

saveRDS(enzymes_ratio, "clean/16_enzymes_ratio.RDS")
