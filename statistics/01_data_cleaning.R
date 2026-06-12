
library(here)

source(here("0_general_functions.R"))

# load data

meta <- read_tsv("data/meta.txt")
dry_matter <- read_tsv("data/dry_matter.txt") %>%
  mutate(sampleid = as.character(sampleid))

# load nmr data
nmr_meta <- read_tsv("data/nmr_meta.txt")
nmr_weighin <- read_tsv("data/nmr_weighin.txt") %>%
  mutate(sampleid = str_remove(sampleid, "AP2_"))
nmr_ileum <- read_tsv("data/nmr_concentrations_ileum.txt", skip = 2)
nmr_feces <- read_tsv("data/nmr_concentrations_feces.txt", skip = 2)

# clean meta

nmr_meta <- nmr_meta %>%
  rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
  dplyr::rename(sampleid = probennummer, nmrid = ...5, ph = ph_probe) %>%
  mutate(sampleid = str_remove(sampleid, "AP2_"))

meta1 <- meta %>%
  mutate(across(c(sampleid, square, animal, period, diet), as.character)) %>%
  left_join(dplyr::select(nmr_meta, sampleid, ph), by = "sampleid") %>% #add ph to meta
  saveRDS("clean/meta1.RDS")

# clean metabolomics data
clean_nmr <- function(input) {
  output <- input %>%
    rename_all(.funs = ~ gsub("\\s+","_", .) %>% tolower) %>%
    dplyr::rename(nmrid = ...1, "dss-d6" = "dss-d6_(chemical_shape_indicator)") %>%
    mutate(nmrid = str_remove(nmrid, "_ex1.cnx")) %>%
    inner_join(dplyr::select(nmr_meta, sampleid, nmrid), by = "nmrid") %>%
    relocate(sampleid) %>%
    dplyr::select(-c(nmrid, "dss-d6")) %>% # remove TSP standard and ethanol
    arrange(sampleid)
  return(output)
}

nmr_ileum_clean <- clean_nmr(nmr_ileum)
saveRDS(nmr_ileum_clean, "clean/5_nmr_ileum_wide_raw.RDS")
nmr_feces_clean <- clean_nmr(nmr_feces)
saveRDS(nmr_feces_clean, "clean/5_nmr_feces_wide_raw.RDS")

# calculate per FM
# original unit = mM = mmol/l = mmol/kg
# account for 0.5 mM concentration of TSP, but scaling to 5 mM -> /10
# account for addition of 60 µl buffer to 54o µl sample -> 
# to calculate per kg FM -> dilution = (weighin+1200µl)/weighin
# mmol/g * dilution
# -> mmol/kg FM

nmr_to_fm <- function(input) {
  output <- input %>%
    pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
    filter(!str_detect(metabolite, "malonate|ethanol")) %>% # exclude (ethanol) and malonate from analysis
    left_join(nmr_weighin, by = "sampleid") %>%
    mutate(concentration = ifelse(is.na(concentration), 0, concentration), # impute missing values with 0
           concentration = concentration/10, # correct for wrong TSP concentration
           concentration = concentration/0.9, # account for buffer addition
           dilution = (1200+weigh_in)/weigh_in,
           concentration = concentration * dilution) %>%
    dplyr::select(-c(dilution, weigh_in)) %>%
    pivot_wider(names_from = "metabolite", values_from = "concentration")
  return(output)
}

nmr_ileum_fm <- nmr_to_fm(nmr_ileum_clean)
saveRDS(nmr_ileum_fm, "clean/5_nmr_ileum_wide_fm.RDS")
write_tsv(nmr_ileum_fm, "tables/5_nmr_ileum_fm.txt")
nmr_feces_fm <- nmr_to_fm(nmr_feces_clean)
saveRDS(nmr_feces_fm, "clean/5_nmr_feces_wide_fm.RDS")
write_tsv(nmr_feces_fm, "tables/5_nmr_feces_fm.txt")

# calculate per dm
# /(DM%/100)
# -> mmol/g DM

nmr_to_dm <- function(input) {
  output <- input %>%
    pivot_longer(-sampleid, names_to = "metabolite", values_to = "concentration") %>%
    left_join(dry_matter, by = "sampleid") %>%
    mutate(concentration = concentration / (dry_matter/100)) %>%
    dplyr::select(-dry_matter) %>%
    pivot_wider(names_from = "metabolite", values_from = "concentration")
  return(output)
}

nmr_ileum_dm <- nmr_to_dm(nmr_ileum_fm)
saveRDS(nmr_ileum_dm, "clean/5_nmr_ileum_wide_dm.RDS")
write_tsv(nmr_ileum_dm, "tables/5_nmr_ileum_dm.txt")
nmr_feces_dm <- nmr_to_dm(nmr_feces_fm)
saveRDS(nmr_feces_dm, "clean/5_nmr_feces_wide_dm.RDS")
write_tsv(nmr_feces_dm, "tables/5_nmr_feces_dm.txt")
