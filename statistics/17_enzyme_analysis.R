
library(here)

source(here("10_enzyme_functions.R"))
source("E:/R/source/ggplot2_theme.R")

# load data

enzymes <- readRDS("clean/16_enzymes_dm.RDS")

meta <- readRDS("clean/meta1.RDS") %>%
  filter(sampleid != "103") # exclude sample 103

# saveswitch
save = FALSE
#save = TRUE

# combine

enz <- meta %>%
  inner_join(enzymes, by = "sampleid")

enz_flow <- rbind(enz_flow_il, enz_flow_fa)

ratio <- meta %>%
  inner_join(enzymes_ratio, by = "sampleid")

# dry_matter

combined_comparison_enzyme(enz, selected_response = "dry_matter", selected_matrix = "faeces",
                    transformation = "tukey", type = "")

combined_comparison_enzyme(enz, selected_response = "dry_matter", selected_matrix = "ileal digesta",
                    transformation = "gauss", type = "hh")

# for cara in faeces

cara_fa <- combined_comparison_and_results(enz, selected_response = "cara_dm", selected_matrix = "faeces", 
                    transformation = "log", type = "hh", response_name = "Carboxypeptidase A", 
                    y_axis = "enzyme activity [U/g DM]", save_name= "cara_fa", save = save)

combined_comparison_enzyme(enz, selected_response = "cara_fm", selected_matrix = "faeces", 
                    transformation = "log", type = "s")

# for cara in ileal digesta

cara_il <- combined_comparison_and_results(enz, selected_response = "cara_dm", selected_matrix = "ileal digesta", 
                    transformation = "tukey", type = "s", response_name = "Carboxypeptidase A", 
                    y_axis = "enzyme activity [U/g DM]", save_name = "cara_il", save = save)

combined_comparison_enzyme(enz, selected_response = "cara_fm", selected_matrix = "ileal digesta", 
                    transformation = "tukey2", type = "s")


# for carb in faeces

carb_fa <- combined_comparison_and_results(enz, selected_response = "carb_dm", selected_matrix = "faeces", 
                    transformation = "log", type = "s", response_name = "Carboxypeptidase B",
                    y_axis = "enzyme activity [U/g DM]", save_name = "carb_fa", save = save)

combined_comparison_enzyme(enz, selected_response = "carb_fm", selected_matrix = "faeces", 
                    transformation = "tukey", type = "s")

# for carb in ileal digesta

carb_il <- combined_comparison_and_results(enz, selected_response = "carb_dm", selected_matrix = "ileal digesta", 
                    transformation = "tukey2", type = "s", response_name = "Carboxypeptidase B", 
                    y_axis = "enzyme activity [U/g DM]", save_name = "carb_il", save = save)

combined_comparison_enzyme(enz, selected_response = "carb_fm", selected_matrix = "ileal digesta", 
                    transformation = "log", type = "s")

# for trypsin in faeces

try_fa <- combined_comparison_and_results(enz, selected_response = "try_dm", selected_matrix = "faeces", 
                    transformation = "tukey2", type = "s", response_name = "Trypsin", 
                    y_axis = "enzyme activity [U/g DM]", save_name = "try_fa", save = save)

combined_comparison_enzyme(enz, selected_response = "try_fm", selected_matrix = "faeces", 
                    transformation = "tukey2", type = "hh")

# for trypsin in ileal digesta

try_il <- combined_comparison_and_results(enz, selected_response = "try_dm", selected_matrix = "ileal digesta", 
                    transformation = "tukey", type = "s", response_name = "Trypsin", 
                    y_axis = "enzyme activity [U/g DM]", save_name = "try_il", save = save)


combined_comparison_enzyme(enz, selected_response = "try_fm", selected_matrix = "ileal digesta", 
                    transformation = "tukey", type = "s")

# for chymotrypsin in faeces

chy_fa <- combined_comparison_and_results(enz, selected_response = "chy_dm", selected_matrix = "faeces", 
                    transformation = "gauss", type = "s", response_name = "Chymotrypsin",
                    y_axis = "enzyme activity [U/g DM]", save_name = "chy_fa", save = save)

combined_comparison_enzyme(enz, selected_response = "chy_fm", selected_matrix = "faeces", 
                    transformation = "tukey2", type = "h")

# for chymotrypsin in ileal digesta

chy_il <- combined_comparison_and_results(enz, selected_response = "chy_dm", selected_matrix = "ileal digesta", 
                    transformation = "none", type = "s", response_name = "Chymotrypsin",
                    y_axis = "enzyme activity [U/g DM]", save_name = "chy_il", save = save)

combined_comparison_enzyme(enz, selected_response = "chy_fm", selected_matrix = "ileal digesta", 
                    transformation = "none", type = "s")

# for amylase in faeces

#amy_fa <- combined_comparison_and_results(enz, selected_response = "amy_dm", selected_matrix = "faeces", 
 #                   transformation = "gauss", type = "s", response_name = "Amylase", 
  #                  y_axis = "enzyme activity [U/g DM]", save_name = "amy_fa", save = save)

script_path <- rstudioapi::getSourceEditorContext()$path
script_name <- basename(script_path)
number <- str_extract(script_name, "[0-9]+")

out <- combined_comparison_enzyme(df = enz, selected_response = "amy_dm", selected_matrix="faeces", 
                                  transformation = "gauss", type = "s")
if (isTRUE(save)) {
  save_big(paste0(number, "_", "amy_fa"))
}

out$cld2$emmean <- out$cld2$emmean*100
out$cld2$SE <- out$cld2$SE*100

amy_fa <- create_results_table(input_df = enz, comparison_object = out, response_name = "Amylase")

create_results_plot(input_df = enz, comparison_object = out, response_name = "Amylase", 
                    y_axis = "enzyme activity [mU/g DM]")
if (isTRUE(save)) {
  save_big(paste0(number, "_result_", "amy_fa"))
}


combined_comparison_enzyme(enz, selected_response = "amy_fm", selected_matrix = "faeces", 
                    transformation = "tukey2", type = "s")

# for amylase in ileal digesta

amy_il <- combined_comparison_and_results(enz, selected_response = "amy_dm", selected_matrix = "ileal digesta", 
                    transformation = "tukey2", type = "hh", response_name = "Amylase", 
                    y_axis = "enzyme activity [U/g DM]", save_name = "amy_il", save = save)

combined_comparison_enzyme(enz, selected_response = "amy_fm", selected_matrix = "ileal digesta", 
                    transformation = "tukey2", type = "hh")

# combine tables

table_il <- cbind(cara_il, carb_il[2], try_il[2], chy_il[2], amy_il[2])
table_fa <- cbind(cara_fa, carb_fa[2], try_fa[2], chy_fa[2], amy_fa[2])
print(table_il)
print(table_fa)
write_tsv(table_il, "tables/17_enzymes_ileum.txt")
write_tsv(table_fa, "tables/17_enzymes_faeces.txt")
save_results_table(table_il, title = "Ileal digesta", name = "17_enzymes_table_il")
save_results_table(table_fa, title = "Faeces", name = "17_enzymes_table_fa")

#save_output()



 
