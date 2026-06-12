rm(list = ls())

setwd("E:/R/R_PiFerm2")

library(tidyverse)

# capture output

#sink(paste0("output/", "00_experimental_design_", str_extract(Sys.time(), "\\d+.\\d+.\\d+"), ".txt"))

######################### randomizing the latin square design

latinsquare = matrix(data = c("diet1", "diet2", "diet3", "diet4", "diet1", "diet2", "diet3", "diet4",
                              "diet4", "diet1", "diet2", "diet3", "diet4", "diet1", "diet2", "diet3",
                              "diet3", "diet4", "diet1", "diet2", "diet3", "diet4", "diet1", "diet2", 
                              "diet2", "diet3", "diet4", "diet1", "diet2", "diet3", "diet4", "diet1"), 
                     nrow = 4, ncol = 8, byrow = T)
latinsquare

set.seed(2024059)
random_col_rep1 <- sample(ncol(latinsquare)/2)
random_col_rep2 <- sample(c(5:8))

random_row_rep1 <- sample(nrow(latinsquare))
random_row_rep2 <- sample(nrow(latinsquare))

latinsquare_random_square1 <- latinsquare[random_row_rep1, random_col_rep1]
latinsquare_random_square1
latinsquare_random_square2 <- latinsquare[random_row_rep2, random_col_rep2]
latinsquare_random_square2

latinsquare_random <- cbind(latinsquare_random_square1, latinsquare_random_square2)
latinsquare_random