
library(here)
library(tidyverse)

source(here("0_general_functions.R"))


# function combining everything for enzymes

combined_comparison_enzyme <- function(df, selected_response, selected_matrix, transformation = "none", 
                                type = c("h", "hh", "s"), model_style = "klein") { #, correction_filter
  filtered_df <- select_response(df, selected_response = selected_response, selected_matrix = selected_matrix) #, correction_filter
  
  out <- combined_comparison(df = filtered_df, transformation = transformation, type = type, model_style = model_style)
  
  plot_pairwise(filtered_df, out, selected_response, selected_matrix, save = save) #, correction_filter
  
  return(out)
}

combined_comparison_and_results <- function(df, selected_response, selected_matrix, transformation = "none",
                                            type = c("h", "hh", "s"), response_name, y_axis,
                                            save_name, save) {
  script_path <- rstudioapi::getSourceEditorContext()$path
  script_name <- basename(script_path)
  number <- str_extract(script_name, "[0-9]+")
  
  out <- combined_comparison_enzyme(df = df, selected_response = selected_response, selected_matrix=selected_matrix, 
                                       transformation = transformation, type = type)
  if (isTRUE(save)) {
    save_big(paste0(number, "_", save_name))
  }
  
  table <- create_results_table(input_df = df, comparison_object = out, response_name = response_name)
  
  create_results_plot(input_df = df, comparison_object = out, response_name = response_name, 
                      y_axis = y_axis)
  if (isTRUE(save)) {
    save_big(paste0(number, "_result_", save_name))
  }
  return(table)
}
