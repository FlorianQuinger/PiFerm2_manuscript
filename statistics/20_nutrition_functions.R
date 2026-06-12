
library(here)
library(tidyverse)

source(here("0_general_functions.R"))



# function combining everything

combined_comparison_nutrition <- function(df, selected_response, transformation = "none", 
                                       type = c("h", "hh", "s")) { #, correction_filter
  filtered_df <- select_response(df, selected_response = selected_response) #, correction_filter
  
  out <- combined_comparison(df = filtered_df, transformation = transformation, type = type)
  
  plot_pairwise(filtered_df, out, selected_response, save = save) #, correction_filter
  
  return(out)
}

# one liner function for each response

combined_comparison_and_results <- function(df, selected_response, transformation = "none",
                                            type = c("h", "hh", "s"), response_name, y_axis,
                                            save_name, save) {
  
  number <- get_script_number()
  
  out <- combined_comparison_nutrition(df = df, selected_response = selected_response, 
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
