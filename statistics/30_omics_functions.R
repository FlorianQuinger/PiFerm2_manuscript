library(vegan)
library(ANCOMBC)
library(patchwork)
library(ROTS)
library(edgeR)
library(MicrobiomeProfiler)
library(clusterProfiler)
library(org.Ss.eg.db)
library(ggrepel)

source(here("0_general_functions.R"))

##################################
# Data cleaning functions
##################################

# sum up to the respective ranks function for data cleaning, to create kreport style table from mmseqs2

calculate_rank_abundance <- function(df, sum_column) {
  ranks = c("R1", "P", "C", "O", "F", "G", "S")
  rank_list <- list()
  for (i in 1:length(ranks)) {
    rank <- ranks[i]
    filtered_df <- df %>%
      dplyr::select(sampleid, !!sym(rank), !!sym(sum_column)) %>%
      filter(!is.na(!!sym(rank))) %>%
      group_by(sampleid, !!sym(rank)) %>%
      summarize(!!sym(sum_column) := sum(!!sym(sum_column)), .groups = "drop") %>%
      pivot_longer(-c(sampleid, !!sym(sum_column)), names_to = "rank", values_to = "name") %>%
      pivot_wider(names_from = sampleid, values_from = !!sym(sum_column)) %>% # introducing NAs that will be replaced by zeros
      pivot_longer(-c(rank, name), names_to = "sampleid", values_to = quo_name(sum_column)) %>%
      mutate(!!sym(sum_column) := ifelse(is.na(!!sym(sum_column)), 0, !!sym(sum_column))) %>% # replace NA with 0
      dplyr::select(name, rank, sampleid, !!sym(sum_column)) %>% # reorder similar to k2report
      {if (sum_column == "rel_abd") group_by(., sampleid) %>% mutate(., rel_abd=rel_abd/sum(rel_abd)*100) %>% ungroup(.) else .} # if rel_abd is choosen, recalculate relative abundances
    
    rank_list[[i]] <- filtered_df
  }
  kstyle_df <- bind_rows(rank_list)
  
  return(kstyle_df)
}

# filter out features with low frequency, takes dataframe with sampleid, abundance, features, and optional rank
# grouping cols should be cols used for frequency calculation -> features and optional rank
# abundance column is rel_abd or reads

filter_by_frequency <- function(df, grouping_cols, abundance_column, frequency_cutoff = 1/3) {
  df_freq <- df %>%
    mutate(abu = ifelse(!!sym(abundance_column) > 0, 1, 0)) %>%
    group_by(!!!syms(grouping_cols)) %>% # syms accepts vector or single string
    mutate(frequency = sum(abu)/length(abu)) %>%
    ungroup() 
  
  p <- ggplot(df_freq, aes(x = frequency)) +
    geom_histogram() +
    #scale_x_log10() +
    geom_vline(xintercept = frequency_cutoff, col = "red", lwd = 1) +
    labs(title = "distribution of frequency values", x = "relative abundance [%/100]")
  print(p)
  
  df_freq_filtered <- df_freq%>%
    filter(frequency > frequency_cutoff) %>%
    dplyr::select(-c(abu, frequency))
  
  print(str_c("Frequency filtering kept", length(unique(df_freq_filtered[[1]])), "of", 
              length(unique(df_freq[[1]])), "features, discarding",
              sum(filter(df_freq, frequency < frequency_cutoff)[[abundance_column]]) 
              / sum(df_freq[[abundance_column]]) *100, 
              "% of", abundance_column,
              sep = " "))
  
  return(df_freq_filtered)
}

# filter out low abundant features, takes dataframe with sampleid, abundance, features, and optional rank
# grouping cols should be cols used for abundance average calculation -> features and optional rank
# abundance column is rel_abd or reads


filter_by_abundance <- function(df, grouping_cols, abundance_column, rank_column = NULL, abundance_cutoff = 0.01) {
  abundance_cutoff <- abundance_cutoff / 100 # convert % into number
  
  df_abu <- df %>%
    {if (is.null(rank_column)) group_by(., sampleid) else group_by(., sampleid, !!sym(rank_column)) } %>%
    mutate(sum = sum(!!sym(abundance_column))) %>%
    ungroup() %>%
    mutate(rel_abd2 = !!sym(abundance_column) / sum) %>%
    group_by(!!!syms(grouping_cols)) %>%
    mutate(avg = mean(rel_abd2)) %>%
    ungroup()
  
  p <- ggplot(df_abu, aes(x = avg)) +
    geom_histogram() +
    scale_x_log10() +
    geom_vline(xintercept = abundance_cutoff, col = "red", lwd = 1) +
    labs(title = "distribution of relative abundance values", x = "relative abundance [%/100]")
  print(p)
  
  df_abu_filtered <- df_abu %>%
    filter(avg > abundance_cutoff) %>% # filtering out features that are in avg lower than the cutoff relative to the sum of abundance at this level
    dplyr::select(-c(sum, rel_abd2, avg))
  
  print(str_c("Abundance filtering kept", length(unique(df_abu_filtered[[1]])), "of", 
              length(unique(df_abu[[1]])), "features, discarding",
              sum(filter(df_abu, avg < abundance_cutoff)$avg) / sum(df_abu$avg) *100, 
              "% of", abundance_column,
              sep = " "))
  
  return(df_abu_filtered)
}

# recalculates rel_abd after some features were filtered out, always calculates to 100 % at each level
# rank column if ktable style is included for grouping

recalculate_rel_abd <- function(df, rank_column = NULL) {
  df_recalc <- df %>%
    {if (is.null(rank_column)) group_by(., sampleid) else group_by(., sampleid, !!sym(rank_column)) } %>%
    mutate(rel_abd = rel_abd / sum(rel_abd) * 100) %>%
    ungroup()
  return(df_recalc)
}

# function combining abundance and frequency filtering for ktable style table and optionally recalculates rel_abd

filter_frequency_and_abundance_ktable <- function(ktable, frequency_cutoff = 1/3, abundance_cutoff = 1e-3) {
  # detect columns: should be features, rank, sampleid, abundance
  feature_column <- colnames(ktable)[1]
  rank_column <- colnames(ktable)[2]
  abundance_column <- colnames(ktable)[4]
  
  print(paste("Feature column is", feature_column))
  print(paste("Rank column is called", rank_column)) 
  print(paste("Abundance column is called", abundance_column))
  
  # process separately for ileum and faeces
  ktable_il <- filter_ileum(ktable)
  ktable_fa <- filter_faeces(ktable)
  
  # filter by frequency
  ktable_il_freq <- filter_by_frequency(df = ktable_il, grouping_cols = c(feature_column, rank_column),
                                        abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
  ktable_fa_freq <- filter_by_frequency(df = ktable_fa, grouping_cols = c(feature_column, rank_column),
                                        abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
  
  #filter out low abundant features
  ktable_il_abu <- filter_by_abundance(df = ktable_il_freq, grouping_cols = c(feature_column, rank_column),
                                       abundance_column = abundance_column, rank_column = rank_column,
                                       abundance_cutoff = abundance_cutoff)
  ktable_fa_abu <- filter_by_abundance(df = ktable_fa_freq, grouping_cols = c(feature_column, rank_column),
                                       abundance_column = abundance_column, rank_column = rank_column,
                                       abundance_cutoff = abundance_cutoff)
  
  # recalculate rel_abd if not reads
  if (abundance_column == "rel_abd") {
    ktable_il_out <- recalculate_rel_abd(df = ktable_il_abu, rank_column = rank_column)
    ktable_fa_out <- recalculate_rel_abd(df = ktable_fa_abu, rank_column = rank_column)
  } else {
    ktable_il_out <- ktable_il_abu
    ktable_fa_out <- ktable_fa_abu
  }
  
  ktable_out <- rbind(ktable_il_out, ktable_fa_out)
  
  return(ktable_out)
}

# combines frequency and abundance filtering for table without ranks

filter_frequency_and_abundance <- function(df, feature_columns = NULL,
                                           frequency_cutoff = 1/3, abundance_cutoff = 1e-3) {
  # detect columns: should be features, sampleid, abundance
  if (is.null(feature_columns)) { # if not manually defined, try to find 
    feature_column <- colnames(df)[1]
    abundance_column <- colnames(df)[3]
  } else {
    feature_column <- feature_columns
    abundance_column <- colnames(df)[length(feature_columns)+2] # number of feature columns determines position of abundance
  }
  
  print(paste("Feature column is", feature_column))
  print(paste("Abundance column is called", abundance_column))
  
  # process separately for ileum and faeces
  df_il <- filter_ileum(df)
  df_fa <- filter_faeces(df)
  
  # filter by frequency
  df_il_freq <- filter_by_frequency(df = df_il, grouping_cols = feature_column,
                                        abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
  df_fa_freq <- filter_by_frequency(df = df_fa, grouping_cols = feature_column,
                                        abundance_column = abundance_column, frequency_cutoff = frequency_cutoff)
  
  #filter out low abundant features
  df_il_abu <- filter_by_abundance(df = df_il_freq, grouping_cols = feature_column,
                                       abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
  df_fa_abu <- filter_by_abundance(df = df_fa_freq, grouping_cols = feature_column,
                                       abundance_column = abundance_column, abundance_cutoff = abundance_cutoff)
  
  # recalculate rel_abd if not reads or intensity
  if (abundance_column == "rel_abd") {
    df_il_out <- recalculate_rel_abd(df_il_abu)
    df_fa_out <- recalculate_rel_abd(df_fa_abu)
  } else {
    df_il_out <- df_il_abu
    df_fa_out <- df_fa_abu
  }
  
  df_out <- rbind(df_il_out, df_fa_out)
  
  return(df_out)
}

# calculate abundances of single ko's from an eggnog style table

calculate_kegg_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), kegg_ko) %>%
    filter(kegg_ko != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(kegg_ko, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(kegg_ko)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, kegg_ko) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    mutate(kegg_ko = str_remove(kegg_ko, "ko:")) %>%
    dplyr::select(kegg_ko, sampleid, !!sym(abundance_column))
  return(output)
}

# calculate abundances of single go's from an eggnog style table, go column should be named "gos"

calculate_go_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), gos) %>%
    filter(gos != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(gos, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(gos)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, gos) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(go = gos, sampleid, !!sym(abundance_column))
  return(output)
}

# caclulate abundance of single cogs from an eggnog style table, cog column should be called cog_accession

calculate_cog_abundance <- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), cog_accession) %>%
    filter(!is.na(cog_accession)) %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(cog_accession, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(cog_accession)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, cog_accession) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(cog_accession, sampleid, !!sym(abundance_column))
  return(output)
}

# calculate abundance of cazymes from eggnog style table

calculate_cazy_abundance<- function(input, abundance_column) {
  output <- input %>%
    dplyr::select(sampleid, !!sym(abundance_column), cazy) %>%
    filter(cazy != "-") %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(cazy, sep = ",") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(cazy)) %>% # calculating factor
    ungroup() %>%
    mutate(!!sym(abundance_column) := !!sym(abundance_column) * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    group_by(sampleid, cazy) %>%
    summarize(!!sym(abundance_column) := sum(!!sym(abundance_column)), .groups = "drop") %>%
    dplyr::select(cazy, sampleid, !!sym(abundance_column))
  return(output)
}

####################################
# data analysis functions
##################################

# function to do pcoa

create_distance_matrix <- function(input_df, region = NULL, dissimilarity_index = "bray", 
                                   abundance_column = "rel_abd") {
  first_column <- colnames(input_df)[1]
  
  if (is.null(region)) {
    print("Warning uncommon columns will be excluded")
  } else if (region == "ileum") {
    input_df <- filter(input_df, str_detect(sampleid, "^1"))
  } else if (region == "faeces") {
    input_df <- filter(input_df, str_detect(sampleid, "^3"))
  }
  
  mat <- input_df %>%
    dplyr::select(!!sym(first_column), sampleid, !!sym(abundance_column)) %>%
    pivot_wider(names_from = first_column, values_from = abundance_column) %>%
    dplyr::select(-sampleid) %>%
    as.matrix()
  
  labels <- input_df %>%
    dplyr::select(!!sym(first_column), sampleid, !!sym(abundance_column)) %>%
    pivot_wider(names_from = first_column, values_from = abundance_column) %>%
    pull(sampleid)
  
  rownames(mat) <- labels
  
  bray <- vegdist(mat, method = dissimilarity_index, na.rm = T)
  
  return(bray)
}

create_pcoa <- function(bray, meta = meta, title = "", second_indicator = "diet") {
  
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  
  pco1 <- paste0("PCo1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100,1), "%)")
  pco2 <- paste0("PCo2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100,1), "%)")
  
  x <- as_tibble(pcoa$points, rownames = "sampleid") %>%
    inner_join(meta, by = "sampleid")
  
  # calculate centroids
  x_centers <- x %>%
    group_by(matrix, diet) %>%
    summarise(c1 = mean(V1), c2 = mean(V2), .groups = "drop")
  
  x_star <- x %>%
    group_by(diet) %>%
    mutate(c1 = mean(V1), c2 = mean(V2)) %>%
    ungroup()
  
  p <- ggplot(x) +
    geom_point(data = x_centers, mapping = aes(x = c1, y = c2,
                                               #color = diet, 
                                               shape = diet, 
                                               fill = matrix),
               show.legend = FALSE, size = 8, alpha = 0.5) +
    #geom_segment(data = x_star, mapping = aes(x = V1, y = V2,
    #                                         xend = c1, yend = c2),
    #           show.legend = F, lwd = 1.5, alpha = 0.2) +
    geom_point(aes(x = V1, y = V2, shape = diet, color = !!sym(second_indicator)), size = 4) +
    #geom_label(aes(label = sampleid)) +
    #geom_label(aes(label = sampling_time)) +
    #geom_label(aes(label = animal)) +
    stat_ellipse(aes(x = V1, y = V2, shape = diet, color = diet, fill = matrix),
                 level = 0.9, lwd = 1.5, alpha = 0.3, show.legend = F) +
    labs(x = pco1, y = pco2, title = title) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(15,16,17,18,21,22,23,24))
  
  return(p)
}

# function to create permutation matrix, needed when observations are unbalanced, expects meta with period column which is used for permuting

create_permutation_matrix <- function(bray, meta = meta, region, treatment_col, design = c("double_latin_square", "row_column")) {
  period_no <- length(unique(meta$period))
  # meta with all sampleids of one region
  if (region == "ileum") {
    filtered_meta <- filter_ileum(meta)
  } else if (region == "faeces") {
    filtered_meta <- filter_faeces(meta)
  }
  # meta with sampleids in bray
  reduced_meta <- as_tibble(as.matrix(bray), rownames = "sampleid") %>%
    dplyr::select(sampleid) %>%
    inner_join(meta, by = "sampleid")
  
  perms <- allPerms(period_no, control = how(observed = T)) # possible permutations if only periods are permuted
  double_perms <- expand.grid(square1 = seq_len(nrow(perms)),square2 = seq_len(nrow(perms)))[-1,] # excluding observed perm
  perms_double <- matrix(ncol = period_no*2, nrow = nrow(double_perms)) # all possible combinations of perms for double latin square
  for (i in 1:nrow(double_perms)) {
    double_perm <- double_perms[i,]
    perm <- c(perms[double_perm[[1]],], perms[double_perm[[2]],]) # concatenate single perms for double perms
    perms_double[i,] <- perm
  }
  
  if (design == "row_column") {
    perms <- perms[-1,] # exclude observed column from perms
  } else if (design == "double_latin_square") {
    perms <- perms_double
  }

  
  permutations <- matrix(ncol = nrow(reduced_meta), nrow = nrow(perms))
  for (i in 1:nrow(perms)) {
    # permute period according to perms
    if (design == "row_column") {
      # permute identical for both squares
      permutation <- c(which(filtered_meta$period == perms[i, 1]),
                       which(filtered_meta$period == perms[i, 2]),
                       which(filtered_meta$period == perms[i, 3]),
                       which(filtered_meta$period == perms[i, 4])) 
    } else if (design == "double_latin_square") {
      # normal order is s1p1, s2p1, s1p2, s2p2, s1p3, s2p3, s1p4, s2p4
      permutation <- c(which(filtered_meta$period == perms[i, 1] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 5] & filtered_meta$square == "2"),
                       which(filtered_meta$period == perms[i, 2] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 6] & filtered_meta$square == "2"),
                       which(filtered_meta$period == perms[i, 3] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 7] & filtered_meta$square == "2"),
                       which(filtered_meta$period == perms[i, 4] & filtered_meta$square == "1"),
                       which(filtered_meta$period == perms[i, 8] & filtered_meta$square == "2"))
    } else {print("no design selected")}
    if (nrow(filtered_meta) != nrow(reduced_meta)) {
      missing_rows <- setdiff(filtered_meta, reduced_meta)
      for (j in nrow(missing_rows)) {
        missing_row <- missing_rows[j,]
        missing_position <- which(filtered_meta$sampleid == missing_row$sampleid) # find position where row was removed
        permutation <- permutation[-missing_position] # remove removed position from permutation
        permuted_position <- which(permutation == missing_position) # detect where missing position was permuted
        missing_label <- filtered_meta[[treatment_col]][missing_position] # detect which label missing position had
        permutation[permuted_position] <- setdiff(which(filtered_meta[[treatment_col]] == missing_label), missing_position)[1] # replace permuted position with identical label from other position
        permutation <- ifelse(permutation > missing_position, permutation - 1, permutation)# substract 1 from each position > than missing position 
      }
    }
    permutations[i,] <- permutation
  }
  
  # for row column design and balanced design, matrix should be the same as
  #ctrl <- with( #main block, never permute between blocks, only within
  #test_meta, how(plots = Plots(strata = animal),  #plot, permutation defined by within argument
   #                           within = Within(type = "free", constant = T))) #constant leads to same within plot permutation for all plots
  
  return(permutations)
}

# do a betadisper test

do_betadisper <- function(bray, meta = meta, permutation_matrix) {
  #filter meta
  x <- as_tibble(as.matrix(bray), rownames = "sampleid") %>%
    dplyr::select(sampleid) %>%
    inner_join(meta, by = "sampleid")
  
  # perform betadisper
  betadis <- betadisper(bray, group= x$diet, type = "centroid")
  permbetadis <- permutest(betadis, permutations = permutation_matrix, pairwise = T)
  
  return(list(betadis,permbetadis))
}

# plot betadispersion as boxplots

plot_betadisper <- function(betadis, meta = meta, title = NULL) {
  # extract permuation p values
  treatment_comparisons <- betadis[[2]]$pairwise$permuted %>%
    as_tibble(rownames = "comparison") %>%
    separate(comparison, into = c("comp1", "comp2"))
  
  # create p matrix
  p_mat <- matrix(nrow = length(unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))),
                  ncol = length(unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))))
  colnames(p_mat) <- unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))
  rownames(p_mat) <- unique(c(treatment_comparisons$comp1, treatment_comparisons$comp2))
  diag(p_mat) <- 1
  
  for (i in 1:nrow(treatment_comparisons)) {
    p_mat[treatment_comparisons$comp1[i], treatment_comparisons$comp2[i]] <- treatment_comparisons$value[i]
    p_mat[treatment_comparisons$comp2[i], treatment_comparisons$comp1[i]] <- treatment_comparisons$value[i]
  }
  
  # get cld
  cld <- generate_cld(p_mat) %>%
    dplyr::rename(diet = treatments) %>%
    mutate(cld = str_remove(treatments_cld, "\\d"))
  
  # get distances
  dispersion <- betadis[[1]]$distances %>%
    as_tibble(rownames = "sampleid") %>%
    left_join(meta, by = "sampleid") %>%
    left_join(cld, by = "diet") %>%
    group_by(diet) %>%
    mutate(mean_dis = mean(value)) %>%
    ungroup()
  
  P <- betadis[[2]]$tab$`Pr(>F)`[1]
  Print <- round_p_value(P)
  
  p <- ggplot(dispersion, aes(x = diet, y = value, fill = diet)) +
    geom_boxplot(width = .5, outliers = F, show.legend = F, alpha = .4) +
    geom_quasirandom(width = .2, size=2, show.legend = F, pch = 21, alpha = .5) +
    {if (P < 0.05) geom_text(aes(y = mean_dis, label = cld))} +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(y = "Distance to centroid", title = title, 
         subtitle = paste("Permuation test: P", Print))
    
  return(p)
}

do_permanova <- function(bray, meta = meta, permutation_matrix) {
  x <- as_tibble(as.matrix(bray), rownames = "sampleid") %>%
    dplyr::select(sampleid) %>%
    inner_join(meta, by = "sampleid")
  # perform PERMANOVA
  
  permanova <- adonis2(bray~diet, data = x, permutations = permutation_matrix)
  
  
  ############ permanova with full formula and without permutation design
  #permanova <- adonis2(bray ~ diet + animal + period, data = x, permutations = 9999, by = "margin")
  
  return(permanova)
}

#generate cld out of matrix and output treatments with cld letters
generate_cld <- function(p_matrix) {
  # insert and absorb algorithm #
  ###############################
  # Get treatment names
  treatments <- rownames(p_matrix)
  
  # find significant differences
  sig <- list()
  counter = 1
  for (i in 1:length(treatments)) { # the treatment to compare with
    for (j in 1:length(treatments)) { # the treatment of interest
      if (p_matrix[i,j] < 0.05) { # if significant different, add to group
        sig[[counter]] <- sort(c(treatments[[i]], treatments[j]))
        counter = counter + 1
      }
    }
  }
  
  # remove identical comparisons
  sig <- unique(sig)
  
  # initialize dataframe
  df <- data.frame(treatments = treatments, column1 = 1)
  
  if (length(sig) > 0) {
    for (i in 1:length(sig)) {
      t1 <- sig[[i]][1]
      t2 <- sig[[i]][2]
      pos_t1 <- which(treatments == t1)
      pos_t2 <- which(treatments == t2)
      # insert
      df_inserted <- df[1]
      for (j in 2:ncol(df)) {
        col = df[,j]
        if (col[pos_t1] == 1 & col[pos_t2] == 1) {
          col1 = col
          col1[pos_t1] <- 0 # replace with 0
          col2 = col
          col2[pos_t2] <- 0 # replace with 0
          df_inserted <- cbind(df_inserted, col1, col2)
        } else {
          df_inserted <- cbind(df_inserted, col) # if no second column needs to be created
        }
      }
      # absorb
      df_absorbed <- df[1]
      for (j in 2:ncol(df_inserted)) {
        contained = 0
        col1 = df_inserted[,j]
        for (k in 2:ncol(df_inserted)) {
          if (j == k) {
            next
          } else {
            col2 = df_inserted[,k]
            # check if col1 is contained in col2
            contained_mini = 1
            for (l in 1:length(col1)) { # iterate over column
              if (col1[l] == 1) { # check only 1 positions
                if (col2[l] == 0) { # if col2 is not 1 at the same position
                  contained_mini = 0 # contained is set FALSE
                }
              }
            }
            # if column is contained, set global contained to TRUE
            if (contained_mini == 1) {
              contained = 1
            }
          }
        }
        if (contained == 0) {
          df_absorbed = cbind(df_absorbed, col1)
        }
      }
      df <- df_absorbed # to start loop again
    }
  }
    
  names(df)[2:ncol(df)] <- sort(letters[1:(ncol(df))-1], decreasing = T)
  
  cld <- df[, c(1,sort(c(2:ncol(df)), decreasing = T))]
  # replace 0 and 1 by letters
  for (i in 2:ncol(cld)) {
    colname = names(cld)[i]
    for (j in 1:nrow(cld)) {
      cld[j,i] <- ifelse(cld[j,i] == 1, colname, "")
    }
  }
  
  # combine treatments with letters
  treatments_cld <- c()
  for (i in 1:nrow(cld)) {
    concat <- str_c(cld[i,],collapse = "")
    treatments_cld[i] <- concat
  }

  return(data.frame(treatments, treatments_cld))
}

pairwise_permanova <- function(bray, meta = meta, treatment = "diet", region, design) {
  treatments <- sort(unique(meta[[which(colnames(meta) == treatment)]]))
  
  #create frame to save values
  p_mat <- matrix(data = 1, nrow = length(treatments), ncol = length(treatments), 
                  dimnames = list(treatments,treatments))
  
  #generate different comparison possibilities
  comparisons <- tibble("comb1" = as.character(), "comb2" = as.character())
  for (i in 1:length(treatments)) {
    comb1 <- treatments[i]
    for (j in 1:length(treatments)) {
      comb2 <- treatments[j]
      # "only upper side of the matrix"...more or less
      if (j > i) {
        comp <- tibble("comb1" = comb1, "comb2" = comb2)
        comparisons <- comparisons %>%
          add_row(comp)
      }
    }
  }
  
  bray_matrix <- as.matrix(bray)
  
  comparisons <- add_column(comparisons, p_value = NA, p_adj = NA)
  
  for (i in 1:nrow(comparisons)) {
    current_comp1 <- comparisons$comb1[i]
    current_comp2 <- comparisons$comb2[i]
    
    # filter relevant sampleids
    if (region == "ileum") {
      meta_filtered <- filter_ileum(meta) %>%
        filter(!!sym(treatment) %in% c(current_comp1, current_comp2))
    } else if (region == "faeces") {
      meta_filtered <- filter_faeces(meta) %>%
        filter(!!sym(treatment) %in% c(current_comp1, current_comp2))
    }
    
    bray_filtered <- as.dist(bray_matrix[rownames(bray_matrix) %in% meta_filtered$sampleid,
                                         colnames(bray_matrix) %in% meta_filtered$sampleid])
    
    permutation_matrix <- create_permutation_matrix(bray = bray_filtered, meta = meta_filtered, region = region,
                                                    treatment_col = treatment, design = design) # not in use right now

    meta_reduced <- as_tibble(as.matrix(bray_filtered), rownames = "sampleid") %>%
      dplyr::select(sampleid) %>%
      inner_join(meta_filtered, by = "sampleid")
    
    permanova <- do_permanova(bray_filtered, meta = meta_reduced, permutation_matrix = 
                                with(meta_reduced, how(plots = Plots(strata = animal), 
                                                        within = Within(type = "free", constant = F))))

    comparisons$p_value[i] <- permanova$`Pr(>F)`[1]
  }
  
  # adjust p value using BH
  #comparisons$p_adj <- p.adjust(comparisons$p_value, method = "BH")
  comparisons$p_adj <- comparisons$p_value # unadjusted due to errors
  
  # fill matrix
  for (i in 1:nrow(comparisons)) {
    current_comp1 <- comparisons$comb1[i]
    current_comp2 <- comparisons$comb2[i]
    p_mat[which(rownames(p_mat) == current_comp1),
          which(colnames(p_mat) == current_comp2)] <- comparisons$p_adj[i]
  }
  
  # mirror matrix
  p_mat[lower.tri(p_mat)] <- t(p_mat)[lower.tri(p_mat)] # transpose to ensure right ordering
  print(p_mat)
  
  # add cld
  cld <- generate_cld(p_mat)
  
  return(cld)
}

do_ordination <- function(input_df, region = NULL, title = "", save = FALSE, save_name = NULL,
                          second_indicator = "diet", dissimilarity_index = "bray", abundance_column = "rel_abd") {
  
  bray <- create_distance_matrix(input_df, region = region, dissimilarity_index = dissimilarity_index,
                                 abundance_column = abundance_column)
  
  if (!is.null(region)){
    permutation_matrix <- create_permutation_matrix(bray = bray, meta = meta, region = region, 
                                                    treatment_col = "diet", design = "double_latin_square")
    
    # perform betadisper
    betadis <- do_betadisper(bray, meta = meta, permutation_matrix = permutation_matrix)
    print(betadis)
    
    p <- plot_betadisper(betadis = betadis, meta = meta, title = title)
    print(p)
    
    if (isTRUE(save)) {
      script_no <- get_script_number()
      save_name2 <- paste0(script_no, "_betadis_", save_name)
      save_big(save_name2)
    }
    
    # perform permanova
    permanova <- do_permanova(bray, meta = meta, permutation_matrix = permutation_matrix)
    print(permanova)
    
    # perform pairwise permanova if significant
    if (permanova$`Pr(>F)`[1] < 0.05) {
      cld <- pairwise_permanova(bray, meta = meta, region = region, design = "double_latin_square")
      # replace diet groups with cld's for plotting
      meta_cld <- meta %>%
        inner_join(cld, by = c("diet"="treatments")) %>%
        dplyr::select(-diet) %>%
        dplyr::rename(diet = treatments_cld)
    } else {meta_cld = meta}
  } else {meta_cld = meta}
  
  p <- create_pcoa(bray, meta = meta_cld, title = title, second_indicator = second_indicator) 
  if (!is.null(region)) {
    Print <- ifelse(permanova$`Pr(>F)`[1] < 0.001, paste("< 0.001"), 
                    paste("=", format(round(permanova$`Pr(>F)`[1],3), nsmall = 3)))
    R2 <- format(round(permanova$R2[1], 3),nsmall = 3)
    p <- p + labs(subtitle = bquote("PERMANOVA:"~R^2 == .(R2)*"," ~ italic(P)~.(Print)))
  }
  return(p)
  
  if (isTRUE(save)) {
    script_no <- get_script_number()
    save_name2 <- paste0(script_no, "_pcoa_", save_name)
    save_big(save_name2)
  }
}

# create taxa barplots from kraken structured table

filter_ktable <- function(ktable, meta, selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                          selected_matrix = c("ileal digesta", "faeces")) {
  filtered_table <- ktable %>%
    inner_join(dplyr::select(meta, sampleid, matrix), by = "sampleid") %>%
    {if (selected_rank != "low") filter(., rank == selected_rank) else .} %>%
    {if (selected_rank != "low") dplyr::select(., -rank) else .} %>%
    filter(matrix == selected_matrix) %>%
    dplyr::select(-matrix)
  return(filtered_table)
}

aggregate_low_abundant_taxa <- function(filtered_table, threshold = 1) {
  aggregated_table <- filtered_table %>%
    group_by(name) %>%
    mutate(mean = mean(rel_abd)) %>%
    ungroup() %>%
    mutate(name = ifelse(mean < threshold, "other", name)) %>%
    group_by(sampleid, name) %>%
    summarize(rel_abd = sum(rel_abd), .groups = "drop")
  return(aggregated_table)
}

create_taxa_barplot <- function(table, meta, title = "") {
  plotting_table <- table %>%
    inner_join(meta, by = "sampleid") %>%
    group_by(diet, name) %>%
    summarise(rel_abd = mean(rel_abd), .groups = "drop") %>%
    mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other"))) # fits other at the end of the legend
  
  ggplot(plotting_table, aes(x = diet, y = rel_abd, fill = name)) +
    geom_bar(stat = "identity", position = "stack", width = 0.8) +
    labs(x = "Diet", y = "Relative abundance (%)", fill = "Taxa", title = title) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0,100.01), expand = c(0,0))
}

# plotting function to create individual barplots for each samples, grouped by a grouping factor
# grouping factor should be a vector c() with a grouping factor for the facetting, and one for the x axis label

create_taxa_barplot_individual <- function(table, meta, title = "", grouping_factors) {
  grouping_factor <- grouping_factors[1]
  x_axis <- grouping_factors[2]
  plotting_table <- table %>%
    inner_join(meta, by = "sampleid") %>%
    mutate(sampleid = str_c(!!sym(grouping_factor), sampleid, sep = "_")) %>%
    mutate(name = factor(name, levels = c(unique(name)[-which(unique(name)=="other")], "other"))) # fits other at the end of the legend
  
  ggplot(plotting_table, aes(x = !!sym(x_axis), y = rel_abd, fill = name)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    labs(y = "Relative abundance [%]", fill = "Taxa", title = title, subtitle = grouping_factor) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0,100.1), expand = c(0,0)) + # 100.1 to not exclude rows due to out of range
    facet_grid(cols= vars(!!sym(grouping_factor))) +
    theme(plot.subtitle = element_text(hjust = 0.5))
}

# if grouping_factors is not defined, a barplot will be created on means per diet, otherwise a bar for each sample will be plotted and ordered according to grouping variables

taxa_barplot_from_ktable <- function(ktable, meta, selected_rank = c("P", "C", "O", "F", "G", "S"), 
                                     selected_matrix = c("ileal digesta", "faeces"), title = "", threshold = 1,
                                     save_name, save = F, grouping_factors = NULL) {
  filtered_table <- filter_ktable(ktable = ktable, meta = meta, selected_rank = selected_rank, 
                                  selected_matrix = selected_matrix)
  
  aggregated_table <- aggregate_low_abundant_taxa(filtered_table = filtered_table, threshold = threshold)
  
  if (is.null(grouping_factors)){
    p <- create_taxa_barplot(aggregated_table, meta = meta, title = title)
  } else {
    p <- create_taxa_barplot_individual(aggregated_table, meta = meta, title = title, 
                                        grouping_factors = grouping_factors)
  }
  
  print(p)
  
  if (isTRUE(save)) {
    save_big(name = save_name)
  }
}

# Differential abundance analysis with ancombc

# columns name, rank, sampleid, abundance expected
prepare_table_for_ancombc <- function(input_df) {
  
  feature_column <- colnames(input_df)[1]
  abundance_column <- colnames(input_df)[3]
  
  matrix <- input_df %>%
    pivot_wider(names_from = feature_column, values_from = abundance_column) %>%
    arrange(sampleid) %>%
    column_to_rownames("sampleid")
  return(matrix)
}

filter_meta <- function(meta, selected_matrix = c("ileal digesta", "faeces")) {
  filtered_meta <- meta %>%
    filter(matrix == selected_matrix)
  return(filtered_meta)
}

# extract taxa from output and change to long format

extract_pairwise_taxa_to_list <- function(output) {
  res_pair_list <- list()
  comparisons <- c("diet2", "diet3", "diet4", "diet3_diet2", "diet4_diet2", "diet4_diet3")
  comparisons_renamed <- c("diet2-diet1", "diet3-diet1", "diet4-diet1", "diet3-diet2", "diet4-diet2", "diet4-diet3")
  
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
      mutate(comparison = comparisons_renamed[i])
    res_pair_list[[i]] <- res_pair_filtered
    names(res_pair_list)[i] <- comparisons_renamed[i]
  }
  return(res_pair_list)
}

# waterfall plots of significant taxa for each comparison

create_waterfall_plot_from_list <- function(res_pair_list, q_threshold = .05, lfc_threshold = 0) {
  res_pair_plot_list <- list()
  for (i in 1:length(res_pair_list)) {
    current_df <- res_pair_list[[i]] %>%
      filter(q < q_threshold) %>%
      mutate(indicator_color = case_when(lfc > lfc_threshold & ss == 1 ~ "#33a02c", # color according to lfc and if passed sensitivity test
                                         lfc > lfc_threshold & ss == 0 ~ "#b2df8a",
                                         lfc < -1*lfc_threshold & ss == 1 ~ "#e31a1c",
                                         lfc < -1*lfc_threshold & ss == 0 ~ "#fb9a99",
                                         .default = "#000000")) %>%
      arrange(lfc) %>%
      mutate(taxon = factor(taxon, levels = taxon))
    
    p <- ggplot(current_df, aes(x = lfc, y = taxon, fill = indicator_color)) +
      geom_bar(stat = "identity", show.legend = T, color = "grey20", width = 0.7) +
      geom_errorbar(aes(xmin = lfc - se, xmax = lfc + se), color = "grey20", width = 0.3) +
      labs(x = expression(Log[2]~fold~change), y = "", title = names(res_pair_list)[i]) + 
      theme(axis.text.y = element_text(size = 8),
            plot.title=element_text(size=22, face="plain", hjust= 0.5, vjust=0)) +
      scale_fill_identity() # use colors in column
    
    res_pair_plot_list[[i]] <- p
  }
  
  p_new <- (res_pair_plot_list[[1]] | res_pair_plot_list[[2]] | res_pair_plot_list[[3]]) /
    (res_pair_plot_list[[4]] | res_pair_plot_list[[5]] | res_pair_plot_list[[6]])
  
  print(p_new)
  
  return(res_pair_plot_list)
}

# volcano plots for each comparison
create_volcano_plot_from_list <- function(res_pair_list, title, q_threshold = .05, lfc_threshold = 0) {
  res_pair_long <- bind_rows(res_pair_list) %>%
    mutate(indicator_color = case_when(q < q_threshold & lfc > lfc_threshold & ss == 1 ~ "#33a02c", # color according to lfc and if passed sensitivity test
                                       q < q_threshold & lfc > lfc_threshold & ss == 0 ~ "#b2df8a",
                                       q < q_threshold & lfc < -1*lfc_threshold & ss == 1 ~ "#e31a1c",
                                       q < q_threshold & lfc < -1*lfc_threshold & ss == 0 ~ "#fb9a99",
                                       .default = "#000000"))
  
  p <- ggplot(res_pair_long, aes(x = lfc, y = -log10(p), color = indicator_color, size = indicator_color)) +
    geom_point(show.legend = F) +
    #geom_hline(yintercept = 0) +
    #geom_hline(yintercept = -log10(q_threshold), linetype = 2) +
    geom_vline(xintercept = 0 + lfc_threshold, linetype = 2, color = "grey30") +
    geom_vline(xintercept = 0 - lfc_threshold, linetype = 2, color = "grey30") +
    scale_color_identity() +
    scale_size_manual(values = c(.2,2,2,2,2)) +
    facet_wrap(~comparison, ncol = 3) +
    labs(x = expression(Log[2]~fold~change), 
         y = expression("Significance ("*-log[10]~italic(P)*-value*")"), title = title)
  
  print(p)
  
  return(p)
}

# combined function to perform ancombc


perform_ancombc_and_plot <- function(ktable, meta, selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                                     selected_matrix = c("ileal digesta", "faeces")) {
  
  filtered_table <- filter_ktable(ktable, meta = meta, selected_rank = selected_rank, 
                                  selected_matrix = selected_matrix)
  
  matrix <- prepare_table_for_ancombc(filtered_table)
  
  filtered_meta <- filter_meta(meta = meta, selected_matrix = selected_matrix) %>%
    column_to_rownames("sampleid")
  
  # perform ancombc
  output <- ancombc2(data = matrix, meta_data = filtered_meta, fix_formula = "diet + animal + period", 
                     #rand_formula = "(1|animal) + (1|period)", # not converging with random effects structure
                     p_adj_method = "holm", group = "diet", n_cl = 8, verbose = T, global = T, pairwise = T,
                     taxa_are_rows = F, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)) # lib_cut = 1000 
  
  res_pair_list <- extract_pairwise_taxa_to_list(output = output)
  
  # build save names
  script_no <- get_script_number()
  save_name <- str_c(script_no, "")
  
  # plot
  waterfall_list <- create_waterfall_plot_from_list(res_pair_list)
  save_big(name = str_c(script_no, "waterfall", str_replace(selected_matrix, " ", "_"), selected_rank, sep = "_"))
  
  create_volcano_plot_from_list(res_pair_list, title = str_c(selected_matrix, selected_rank, sep = " "))
  save_big(name = str_c(script_no, "volcano", str_replace(selected_matrix, " ", "_"), selected_rank, sep = "_"))
  
  # create output object
  out <- list()
  out[["input"]] <- filtered_table
  out[["meta"]] <- filtered_meta
  out[["ancombc_output"]] <- output
  out[["res_pair_list"]] <- res_pair_list
  out[["waterfall_list"]] <- waterfall_list
  
  return(out)
}


# helper function to create pval matrix from dataframe with p vals

create_p_matrix_from_df <- function(taxon_pvals) {
  levels <- levels(meta$diet)
  p_matrix <- matrix(nrow = 4, ncol = 4, dimnames = list(levels,
                                                         levels))
  diag(p_matrix) <- 1
  
  for (i in 1:nrow(taxon_pvals)) {
    row <- taxon_pvals[i,] %>%
      separate(comparison, into = c("treat1", "treat2"), sep = "\\s+vs\\.\\s+|-")
    p_matrix[which(rownames(p_matrix) == row$treat1), which(colnames(p_matrix) == row$treat2)] <- row$q
    p_matrix[which(rownames(p_matrix) == row$treat2), which(colnames(p_matrix) == row$treat1)] <- row$q
  }
  
  return(p_matrix)
}

# create filtered plot from ancombc output

create_plot_for_taxon <- function(output_object, df_rel_abd, meta, 
                                  selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                                  selected_matrix = c("ileal digesta", "faeces"), selected_taxon,
                                  save = F) {
  
  filtered_df <- filter_ktable(df_rel_abd, meta = meta, selected_rank = selected_rank, 
                               selected_matrix = selected_matrix)
  
  p_adj <- filter(output_object$ancombc_output$res_global, taxon == selected_taxon)$q_val
  Print <- round_p_value(p_adj)
  
  taxon_pvals <- bind_rows(output_object$res_pair_list) %>%
    filter(taxon == selected_taxon)
  
  p_matrix <- create_p_matrix_from_df(taxon_pvals)
  
  cld <- generate_cld(p_matrix) %>%
    mutate(cld = treatments_cld,
           cld = str_remove(cld, treatments),
           diet = str_remove(treatments, "diet")) 
  
  taxon_df <- filtered_df %>%
    filter(name == selected_taxon) %>%
    left_join(dplyr::select(meta, sampleid, diet), by = "sampleid") %>%
    left_join(cld, by = "diet") %>%
    group_by(diet) %>%
    mutate(mean = mean(rel_abd)) %>%
    ungroup()
  
  p <- ggplot(taxon_df, aes(x = diet, y = rel_abd)) +
    geom_boxplot(outliers = F, fill = "snow2", width = 0.3) +
    geom_beeswarm() +
    geom_text(aes(y = mean, label = cld), size = 5, position = position_nudge(x = 0.2)) +
    annotate("text", x = 2.5, y = max(taxon_df$rel_abd), label = paste("q", Print), size = 6) +
    labs(x = "Diet", y = "Relative abundance (%)", title = selected_taxon)
  print(p)
  
  if (isTRUE(save)) {
    # build save names
    script_no <- get_script_number()
    
    save_big(name = str_c(script_no, "taxon_boxplot", str_replace(selected_matrix, " ", "_"), selected_rank, 
                          selected_taxon, sep = "_"))
  }
}

# loop taxa plot function ofer a vector of significant taxa

create_plots_for_taxa <- function(output_object, df_rel_abd, meta, 
                                  selected_rank = c("P", "C", "O", "F", "G", "S", "low"), 
                                  selected_matrix = c("ileal digesta", "faeces"), 
                                  taxa_vector, save) {
  for (i in 1:length(taxa_vector)) {
    create_plot_for_taxon(output_object = output_object, df_rel_abd = df_rel_abd, meta = meta,
                          selected_rank = selected_rank, selected_matrix = selected_matrix,
                          selected_taxon = taxa_vector[i], save = save)
  }
}


# edgeR

# creates edgeR object from table with sampleid column, feature column in first place, and custom intensity or reads column

create_edger_object <- function(input_df, abundance_column) {
  feature_column <- colnames(input_df)[1]
  
  output <- input_df %>%
    dplyr::select(!!sym(feature_column), sampleid, !!sym(abundance_column)) %>%
    pivot_wider(names_from = sampleid, values_from = !!sym(abundance_column))
  mat <- as.matrix(output[,-1])
  rownames(mat) <- output[[1]]
  
  groups <- tibble("sampleid" = colnames(mat)) %>%
    inner_join(meta, by = "sampleid") %>%
    dplyr::select(sampleid, animal, period, diet)
  # combine to edger object
  edger_object <- DGEList(counts = mat, samples=groups)
  
  return(edger_object)
}

edger_analysis <- function(edger_object, feature_column) {
  #create design matrix
  design <- model.matrix(~diet+animal+period, data = edger_object$samples)
  rownames(design) <- edger_object$samples$sampleid
  #estimate dispersion
  disp <- estimateDisp(edger_object, design)
  plotBCV(disp)
  #fit glm
  fit <- glmFit(disp, design = design)
  
  # do anova
  anova <- glmLRT(fit, coef = 2:4)
  anova_adj <- decideTests(anova)
  anova_adj_sum <- summary(anova_adj)
  
  # continue if significant
  if (anova_adj_sum[2] == 0) {
    print("No significant differences.")
    plotMD(anova)
  } else {
    contrasts <- design[0,] 
    contrasts <- rbind(contrasts,
                       "diet2-diet1" = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0),
                       "diet3-diet1" = c(0,0,1,0,0,0,0,0,0,0,0,0,0,0),
                       "diet4-diet1" = c(0,0,0,1,0,0,0,0,0,0,0,0,0,0),
                       "diet3-diet2" = c(0,-1,1,0,0,0,0,0,0,0,0,0,0,0),
                       "diet4-diet2" = c(0,-1,0,1,0,0,0,0,0,0,0,0,0,0),
                       "diet4-diet3" = c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0))
    results <- vector("list", length = nrow(contrasts))
    names(results) <- rownames(contrasts)
    
    for (i in 1:nrow(contrasts)) {
      contrast <- as.vector(contrasts[i,])
      result <- glmTreat(fit, contrast = contrast)#, lfc = 1) # test + filters lfc > 1
      sig <- decideTests(result, adjust.method = "BH") # adjust P and filter according to 0.05
      table <- topTags(result, n=nrow(fit$counts), adjust.method = "BH", sort.by = "none", p.value = 1)$table # output all unsorted
      print(summary(sig))
      #save in output
      results[[i]][["result"]] <- result
      results[[i]][["significant"]] <- sig
      results[[i]][["table"]] <- table
    }

  }
  return(results)
}

# takes edger output list and filters every table for significant results
# stores output table in edger output object

create_table_sig_from_edger <- function(edger_results, feature_column) {
  for (i in 1:length(edger_results)) {
    table <- edger_results[[i]][["table"]] %>%
      as_tibble(rownames = feature_column) %>%
      mutate(contrast = names(edger_results)[i],
             sig = ifelse(FDR < 0.05, "sig", "not sig"),
             sig = ifelse(is.na(sig), "not sig", sig))
    edger_results[[i]][["table_sig"]] <- table # safe table to list for later
  }
  
  return(edger_results)
}

# concatenates edger outputs and creates a volcano plot

create_volcano_plot_from_edger <- function(edger_results, title = "", no_label = 0) {
  
  # concatenate all results and create plot
  for (i in 1:length(edger_results)) {
    table <- edger_results[[i]][["table_sig"]]
    if (i == 1) {
      all_results <- table 
    } else { # from second iteration onwards combine previous tables with latest table
      all_results <- rbind(all_results, table)
    }
  }
  
  first_column <- colnames(all_results)[1]
  # add labels for ggrepel
  all_results <- all_results %>%
    group_by(contrast, sig) %>%
    mutate(rank_asc = row_number(logFC),
           rank_desc = row_number(-logFC),
           label = ifelse(sig == "sig" & ((rank_desc <= no_label & logFC > 0 )| (rank_asc <= no_label & logFC < 0)), 
                          !!sym(first_column), "")) %>%
    ungroup()
  
  p <- ggplot(all_results, aes(x = logFC, y = -log10(PValue), color = sig, size = sig)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey30") +
    geom_point(show.legend = T) +
    facet_wrap(~contrast, ncol = 3) +
    scale_color_manual(values = c("grey20", "red"), labels = c("q ≥ 0.05", "q < 0.05")) +
    scale_size_manual(values = c(0.1,1), labels = c("q ≥ 0.05", "q < 0.05")) +
    labs(x = expression(Log[2]~fold~change), 
         y = expression("Significance ("*-log[10]~italic(P)*-value*")"), title = title, 
         color = "Significance", size = "Significance") +
    theme(legend.position = "bottom")
    
  if (no_label != 0) {
    p <- p +
      geom_text_repel(aes(x = logFC, y = -log10(PValue), label = label), color = "black", show.legend = F, 
                      max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, max.time = 3, size = 5, 
                      inherit.aes = F, segment.color = "grey30") +
      scale_x_continuous(expand = expansion(mult = 0.2)) +
      scale_y_continuous(expand = expansion(mult = c(0,0.2)))
  }
  return(p)
}

perform_edger_and_plot <- function(input_df, selected_matrix = c("ileum", "faeces"), abundance_column, save_name, save = F) {
  
  # filter
  if (selected_matrix == "ileum") {
    filtered_df <- filter_ileum(input_df)
  } else if (selected_matrix == "faeces") {
    filtered_df <- filter_faeces(input_df)
  }
  
  edger_object <- create_edger_object(filtered_df, abundance_column = abundance_column)
  
  edger_results <- edger_analysis(edger_object = edger_object, 
                                  feature_column = colnames(input_df)[1])
  
  edger_results <- create_table_sig_from_edger(edger_results = edger_results, feature_column = colnames(input_df)[1])
  
  p <- create_volcano_plot_from_edger(edger_results = edger_results, title = str_c(save_name, selected_matrix, sep = " "))
  print(p)
  
  # build save names
  if (isTRUE(save)) {
    script_no <- get_script_number()
    save_big(name = str_c(script_no, save_name, str_replace(selected_matrix, " ", "_"), sep = "_"))
  }
  
  return(edger_results)
}

# add functions to edger result object based on meta file with all functions
# join_by defines column to join, must be same in results object and functions_object

add_functions <- function(results_object, functions_object, join_by) {
  for (i in 1:length(results_object)) {
    table_func <- results_object[[i]][["table_sig"]] %>%
      #filter(sig == "sig") %>%
      left_join(functions_object, by = join_by)
    results_object[[i]][["table_sig"]] <- table_func # save as table sig
  }
  return(results_object)
}

# add taxonomy to edger results object based on meta file with all taxa assignments
# join_by defines join column, but should be bins

add_taxonomy <- function(results_object, taxonomy_object, join_by) {
  for (i in 1:length(results_object)) {
    table_func <- results_object[[i]][["table_sig"]] %>%
      mutate(bin = proteinid) %>%
      separate(bin, into = "bin", sep = "_") %>%
      #filter(sig == "sig") %>%
      left_join(taxonomy_object, by = join_by)
    results_object[[i]][["table_sig"]] <- table_func # save as table sig
  }
  return(results_object)
}

# from protein or gene differential abundance analysis with edger
# relocates the function column to be used for enrichment analysis

transform_to_function <- function(results_object, function_name) {
  for (i in 1:length(results_object)) {
    table_func <- results_object[[i]][["table_sig"]] %>%
      filter(!is.na(!!sym(function_name))) %>%
      filter(!!sym(function_name) != "-") %>%
      mutate(!!sym(function_name) := str_remove(!!sym(function_name), ";$")) %>%
      separate_rows(!!sym(function_name), sep = ",") %>%
      separate_rows(!!sym(function_name), sep = ";") %>%    
      mutate(!!sym(function_name) := str_remove_all(!!sym(function_name), "^ssc:")) %>%
      mutate(!!sym(function_name) := str_remove_all(!!sym(function_name), "^ko:")) %>%
      mutate(!!sym(function_name) := str_remove_all(!!sym(function_name), "\\s")) %>%
      dplyr::relocate(!!sym(function_name))
    
    results_object[[i]][["table_sig"]] <- table_func
  }
  return(results_object)
}

# Enrichment analysis from edger results object
# input should be an edger results object with kegg ko's in first column of table_sig 

enrich_kegg <- function(results_object, save_name, save) {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    feature_name <- names(table_sig)[1] # kegg ko's must be the first column
    universe <- table_sig[[feature_name]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(feature_name))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(feature_name))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichKO(overexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(name = paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(paste("Overexpressed KEGGs for", contrast))
        cat(overexpressed)
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichKO(underexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(paste("Underexpressed KEGGs for", contrast))
        cat(underexpressed)
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# enrichment function for kegg modules
# still a bit buggy
enrich_mkegg <- function(results_object, save_name, save) { 
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    feature_name <- names(table_sig)[1] # kegg ko's must be the first column
    universe <- table_sig[[feature_name]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(feature_name))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(feature_name))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichModule(overexpressed, universe = universe, minGSSize = 6) # minGsize not recognized correct
      if (!is.null(enrichment)) {
        if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
          p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
          print(p)
          if (isTRUE(save)) {
            save_big(name = paste0(save_name, "_", contrast, "_overexpressed"))
          }
          print(paste("Overexpressed KEGGs for", contrast))
          cat(overexpressed)
        } else {
          print(paste("No enriched KEGGs found for", contrast))
        }
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichModule(underexpressed, universe = universe)
      if (!is.null(enrichment)) {
        if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
          p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
          print(p)
          if (isTRUE(save)) {
            save_big(paste0(save_name, "_", contrast, "_underexpressed"))
          }
          print(paste("Underexpressed KEGGs for", contrast))
          cat(underexpressed)
        } else {
          print(paste("No enriched KEGGs found for", contrast))
        }
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# creates a heatmap from an edger object and a object containing only significant proteins or genes incl. rel_abd per sample
# heatmap contains cld

create_heatmap_from_edger <- function(results_object, sig_df, function_df = NULL, heatmap_y = "protein_names",
                                      save_name, save) {
  
  first_column <- colnames(sig_df)[1]
  
  sig_proteins <- unique(sig_df[[sym(first_column)]])
  
  #aggregate results in one table
  for (i in 1:length(results_object)) {
    protein_pval <- results_object[[i]]$table_sig %>%
      dplyr::select(first_column, comparison = contrast, q = FDR) # adapt to previous functions
    if (i == 1) {
      protein_pvals <- protein_pval
    } else {
      protein_pvals <- rbind(protein_pvals, protein_pval)
    }
  }
  
  protein_pvals <- protein_pvals %>%
    filter(!!sym(first_column) %in% sig_proteins)
  
  for (i in 1:length(sig_proteins)) {
    protein_pvals_filtered <- protein_pvals %>%
      filter(!!sym(first_column) == sig_proteins[i])
    
    p_matrix <- create_p_matrix_from_df(protein_pvals_filtered)
    
    cld <- generate_cld(p_matrix) %>%
      mutate(cld = treatments_cld,
             cld = str_remove(cld, treatments),
             diet = str_remove(treatments, "diet")) %>%
      add_column(!!sym(first_column) := sig_proteins[i]) %>%
      dplyr::select(first_column, diet, cld)
    
    if (i == 1) {
      protein_cld <- cld
    } else {
      protein_cld <- rbind(protein_cld, cld)
    }
  }
  
  heatmap <- sig_df %>%
    left_join(meta, by = "sampleid") %>%
    group_by(!!sym(first_column), diet) %>%  # identical protein-names will be averaged
    summarise(rel_abd = median(rel_abd), .groups = "drop") %>% # replaced mean by median
    {if (!is.null(function_df)) left_join(., function_df, by = first_column) else .} %>% 
    left_join(protein_cld, by = c(first_column, "diet")) %>%
    group_by(!!sym(first_column)) %>%
    mutate(value = scale(rel_abd)) %>%
    ungroup()
  
  p <- ggplot(heatmap, aes(x = diet, y = !!sym(heatmap_y), fill = value)) +
    geom_tile(show.legend = F) +
    geom_text(aes(label = cld), size = 2) +
    coord_fixed() +
    scale_y_discrete(position = "right") +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 8)) +
    labs(y = "")
  print(p)
  
  if (isTRUE(save)) {
    save_big(name = save_name)
  }
}


##############################
# Functions from Metaproteomics
##############################

# calculate taxa from protein intensities
# takes protein df and genomes df

calculate_taxa_intensity <- function(input, genomes) {
  output <- input %>%
    mutate(id = row_number()) %>% # factor used for intensity splitting, for every row
    separate_rows(proteingroup, sep = ";") %>%
    group_by(id) %>%
    mutate(factor = 1 / length(proteingroup)) %>% # calculating factor
    ungroup() %>%
    mutate(intensity = intensity * factor) %>% # apply factor to intensities
    dplyr::select(-c(factor, id)) %>%
    filter(str_detect(proteingroup, "^MGYG")) %>% # filter out host proteins
    separate(proteingroup, into = "bin", sep = "_") %>%
    group_by(sampleid, bin) %>%
    summarise(intensity = sum(intensity), .groups = "drop") %>% # sum up intensities per sample and genome
    left_join(genomes, by = "bin") %>%
    dplyr::select(bin, sampleid, intensity, "R1", P, C, O, "F", G, S)
  return(output)
}


# enrichment analysis of COGs for microbial proteins
# takes edger results object
# first column of sig df should be COGs -> either from cog dea analysis or translated to cog

enrich_cog <- function(results_object, save_name, save) {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    first_column <- colnames(table_sig)[1]
    universe <- table_sig[["first_column"]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(first_column))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(first_column))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichCOG(overexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(overexpressed)
      } else {
        print(paste("No enriched COGs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed COGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichCOG(underexpressed, universe = universe)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(underexpressed)
      } else {
        print(paste("No enriched COGs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed COGs for", contrast))
    }
  }
}

# Enrichtment analysis for host keggs
# takes edger results object
# first column of sig df should be keggs -> either from cog dea analysis or translated to kegg
# default organism is sus scrofa, "psat" for pisum sativum

enrich_kegg_host <- function(results_object, save_name, save, organism = "ssc") {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    first_column <- colnames(table_sig)[1]
    universe <- table_sig[[first_column]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(first_column))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(first_column))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichKEGG(overexpressed, universe = universe, organism = organism)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(overexpressed)
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichKEGG(underexpressed, universe = universe, organism = organism)
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(underexpressed)
      } else {
        print(paste("No enriched KEGGs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# Enrichtment analysis for host gos
# takes edger results object
# first column of sig df should be go -> either from cog dea analysis or translated to go

enrich_go_host <- function(results_object, save_name, save) {
  for (i in 1:length(results_object)) {
    contrast <- names(results_object)[i]
    table_sig <- results_object[[i]][["table_sig"]]
    first_column <- colnames(table_sig)[1]
    universe <- table_sig[[first_column]]
    overexpressed <- table_sig %>%
      filter(sig == "sig" & logFC > 0) %>%
      pull(!!sym(first_column))
    underexpressed <- table_sig %>%
      filter(sig == "sig" & logFC < 0) %>%
      pull(!!sym(first_column))
    
    if (length(overexpressed)>=10) {
      enrichment <- enrichGO(overexpressed, "org.Ss.eg.db",universe = universe, keyType = "GO", ont = "BP")
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "overexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_overexpressed"))
        }
        print(overexpressed)
      } else {
        print(paste("No enriched GOs found for", contrast))
      }
    } else {
      print(paste("Not enough overexpressed KEGGs for", contrast))
    }
    
    if (length(underexpressed)>=10) {
      enrichment <- enrichGO(underexpressed, "org.Ss.eg.db",universe = universe, keyType = "GO", ont = "BP")
      if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
        p <- enrichplot::dotplot(enrichment) + labs(title = paste(contrast, "underexpressed"))
        print(p)
        if (isTRUE(save)) {
          save_big(paste0(save_name, "_", contrast, "_underexpressed"))
        }
        print(underexpressed)
      } else {
        print(paste("No enriched GOs found for", contrast))
      }
    } else {
      print(paste("Not enough underexpressed KEGGs for", contrast))
    }
  }
}

# output significant proteins etc for correlation purposes

get_significant <- function(venn_input, input_df) {
  first_column <- colnames(input_df)[1]
  output <- venn_input %>%
    as_tibble(rownames = first_column) %>%
    dplyr::select(!!sym(first_column)) %>%
    inner_join(input_df, by = first_column) 
  return(output)
}

# Venn diagramms
# takes results object of ancombc or edgeR and input_df with rel_abd
# selected_rank and selected_matrix need to be specified in case of a ktable
# selected_matrix must be specified to filter for matrix

create_venn_for_sig <- function(results_object, input_df, selected_matrix = c("ileal digesta", "faeces"),
                                selected_rank = NULL, save_name, save) {
  # pre test for ancombc objects
  if (names(results_object)[1] == "input") { # indicator for ancombc object
    results_object <- results_object[["res_pair_list"]] # change results object to res_pair_list
  }
  
  for (i in 1:length(results_object)) {
    first_column <- colnames(input_df)[1]
    contrast <- names(results_object)[i]
    table <- results_object[[i]]
    #ancombc vs edgeR block
    if (is.data.frame(table)) { # detect if table is already dataframe (ancombc) or list (edgeR)
      table <- table %>%
        mutate(sig = ifelse(q < 0.05, "sig", NA)) %>% # if ancombc add sig column
        dplyr::rename(!!sym(first_column) := taxon) # unify name of first column
    } else {
      table <- table[["table_sig"]] # choose table_sig if edgeR
    }
     table_sig <- table %>%
      dplyr::select(!!sym(first_column), sig) %>%
      mutate(sig = ifelse(sig == "sig", 1, 0),
             sig = ifelse(is.na(sig), 0, sig)) %>% # replace NAs separately
      dplyr::rename(!!sym(contrast):=sig)
    if (i == 1) {
      all_sig <- table_sig
    } else { # from second iteration onwards combine previous tables with latest table
      all_sig <- all_sig %>%
        left_join(table_sig, by = first_column)
    }
  }
  sig_df <- as.data.frame(all_sig[,-1])
  rownames(sig_df) <- all_sig[[first_column]]
  sig_df <- sig_df[rowSums(sig_df) > 0,]
  
  #create venn
  save_venn(sig_df)
  
  # create sig output table
  if (is.null(selected_rank)) { # if no rank is specified -> no ktable
    output_df <- filter_ktable(input_df, meta = meta, selected_rank = "low", 
                              selected_matrix = selected_matrix)
  } else {
    output_df <- filter_ktable(input_df, meta = meta, selected_rank = selected_rank, 
                              selected_matrix = selected_matrix)
  }
  
  output <- get_significant(sig_df, output_df) # return tibble with rel_abd for correlating
  
  # create save name and save
  if (isTRUE(save)) {
    no <- get_script_number()
    matrix <- str_extract(selected_matrix, ".{2}")
    saveRDS(output, file = str_c("clean/", str_c(no, "sig", matrix, selected_rank, save_name, sep = "_"), ".RDS"))
  }
  
  return(output)
}

##############################
# functions from metabolomics
#############################

# function for preparation

prepare_metabolomics <- function(input) {
  output <- input %>%
    inner_join(dplyr::select(meta, sampleid, square, animal, period, diet)) %>%
    relocate(animal, period, square, diet) %>%
    dplyr::select(-sampleid)
  return(output)
}

# analysis function

loop_comparison_metabolimics <- function(input_df, y_axis, table_title, save_name, save = save) {
  shapiro_table <- data.frame("metabolite" = as.character(), "p_value"=as.numeric(), transformation = as.character()) # to check for bad transformations
  pb <- txtProgressBar(min = 1, max = ncol(input_df), style = 3) # initialize progessbar
  for (i in 5:ncol(input_df)) { #assuming 4 meta columns in the df
    setTxtProgressBar(pb, i) # print progessbar
    
    number <- get_script_number()
    
    current_response <- colnames(input_df)[i]
    # filter df for i'th response
    filtered_df <- select_response(input_df, current_response)
    # remove rows with NA (for functional activity)
    filtered_df <- filter(filtered_df, !is.na(response))
    # perform reml
    out <- combined_comparison(df = filtered_df, transformation = "test")
    # save shapiro results
    shapiro_table <- rbind(shapiro_table, c(current_response, round(out$shapiro$p.value, 3), out$transformation[1]))
    # print plot
    plot_pairwise(filtered_df = filtered_df, output = out, selected_response = current_response)
    if (isTRUE(save)) {
      save_big(paste0(number, "_", save_name, "_", current_response))
    }
    # print and save plot
    create_results_plot(input_df = filtered_df, comparison_object = out, response_name = current_response, y_axis = y_axis)
    if (isTRUE(save)) {
      save_big(paste0(number, "_result_", save_name, "_", current_response))
    }
    table <- create_results_table(input_df = filtered_df, comparison_object = out, response_name = current_response, digits = 1)
    table <- rbind(table, c("exact_p", out$ano[which(rownames(out$ano) == "diet"), 
                                                             which(colnames(out$ano) == "Pr(>F)")])) # for exact p value adjustment
    
    if (i == 5) { # create table if first iteration, join afterwards
      final_table <- table
    } else {
      final_table <- final_table %>%
        left_join(table, by = "diet")
    }
    
  }
  final_table <- rbind(final_table[-7,], 
                       c("P-adj", 
                         round(p.adjust(final_table[7, 2:ncol(final_table)], method = "BH"), 3)))
  if (isTRUE(save)) {
    final_table_gt <- save_results_table(final_table, title = table_title, name =  paste0(number, "_table_", save_name))
    gtsave(final_table_gt, filename = paste0(number, "_table_", save_name, ".png"), "plots", vwidth = 3000)
  }
  
  close(pb) # close progessbar
  print(shapiro_table)
  return(final_table)
}


###############################
# Functions for functional activity
#######################################

# function to join a gene and a protein frame by sampleid and join_column2
# df's should be filtered for taxa level beforehand
# gets filtered per matrix for missingness < than missingness and imputes 0s by a small value

create_integrated_df <- function(df1, df2, join_column2, missingness = 0) {
  integrated <- inner_join(df1, df2, by = c("sampleid", join_column2))
  
  # filtering for ileum
  integrated_il <- filter_ileum(integrated) %>%
    group_by(!!sym(join_column2)) %>%
    filter((sum(rel_abd_prot == 0)/length(rel_abd_prot)) <= missingness & 
             (sum(rel_abd_gene == 0)/length(rel_abd_gene)) <= missingness) %>%
    ungroup()
  
  # filtering for faeces
  integrated_fa <- filter_faeces(integrated) %>%
    group_by(!!sym(join_column2)) %>%
    filter((sum(rel_abd_prot == 0)/length(rel_abd_prot)) <= missingness & 
             (sum(rel_abd_gene == 0)/length(rel_abd_gene)) <= missingness) %>%
    ungroup()
  
  # no imputation but removal of Inf
  integrated <- rbind(integrated_il, integrated_fa) %>%
    # imputation with a small value
    # mutate(rel_abd_prot = ifelse(rel_abd_prot == 0, 1e-8, rel_abd_prot), # impute with a low value to mitigate Inf
    #        rel_abd_gene = ifelse(rel_abd_gene == 0, 1e-8, rel_abd_gene)) %>%
    # imputation with minimum value
    # mutate(rel_abd_prot = ifelse(rel_abd_prot == 0, NA, rel_abd_prot), # replace 0 with NA
    #        rel_abd_gene = ifelse(rel_abd_gene == 0, NA, rel_abd_gene)) %>%
    # group_by(!!sym(join_column2)) %>%
    # mutate(rel_abd_prot = ifelse(is.na(rel_abd_prot), min(rel_abd_prot, na.rm = T), rel_abd_prot), # impute minimum without NA
    #        rel_abd_gene = ifelse(is.na(rel_abd_gene), min(rel_abd_gene, na.rm = T), rel_abd_gene)) %>%
    # ungroup() %>%
    mutate(func_act = log2(rel_abd_prot/rel_abd_gene)) %>%
    filter(func_act != Inf) %>%
    filter(func_act != -Inf)
  
  return(integrated)
}

# function for filtering and plotting

filter_results_table <- function(results_table, p_threshold, feature_name) {
  filtered_table <- results_table %>%
    column_to_rownames("diet") %>%
    t() %>%
    as_tibble(rownames = feature_name) %>%
    mutate(`P-adj` = as.numeric(`P-adj`),
           `P-value` = as.numeric(`P-value`)) %>%
    filter(`P-adj` < p_threshold) %>%
    pivot_longer(c("1", "2", "3", "4"), names_to = "diet", values_to = "letter") %>%
    mutate(mean = as.numeric(str_remove_all(letter, "[a-z]+$")),
           letter = str_extract_all(letter, "[a-z]+$"))
  return(filtered_table)
}

plot_functional_activity <- function(integrated_df, sig_results, region = c("ileum", "faeces"), name_column = NULL) {
  feature_column <- colnames(sig_results)[1]
  
  # if name_column not defined, define it as feature column
  if (is.null(name_column)) {
    name_column <- feature_column
  }
  
  if (region == "ileum") {
    filtered_df <- filter_ileum(integrated_df)
  } else {
    filtered_df <- filter_faeces(integrated_df)
  }
  
  plotting_df <- filtered_df %>%
    inner_join(dplyr::select(meta, sampleid, square, animal, period, diet), by = "sampleid") %>%
    inner_join(sig_results, by = c(feature_column, "diet"))
  
  p <- ggplot(plotting_df, aes(x = diet)) +
    facet_wrap(facets = vars(!!sym(name_column)), scales = "free_y") +
    geom_boxplot(aes(y = func_act), fill = "#A4BED5FF", outliers = F, alpha = .5) +
    geom_quasirandom(aes(y = func_act, size = rel_abd_prot), color = "#023743FF", width = 0.2, alpha = .7) +
    geom_quasirandom(aes(y = func_act, size = rel_abd_gene), color = "#FED789FF", show.legend = F, width = 0.2, 
                     alpha = .5) +
    labs(size = "rel. abd. %", y = "log2 protein/gene ratio") +
    geom_text(aes(x = 2.5, y = max(func_act), label = paste("q = ", `P-adj`))) +
    geom_text(aes(y = mean, label = letter))
  print(p)
  return(dplyr::select(plotting_df, !!sym(feature_column), sampleid, func_act))
}

# function taking an integrated frame with columns for rel_abd_gene and rel_abd_prot and precalculated func_act
# frames should be filtered for region
# testing dna agains protein rel_abd with wilcoxon test and creating a volcano plot based on mean func_act and adj p

test_protein_against_dna <- function(input) {
  # remove Inf values
  input <- input %>%
    filter(func_act != Inf)
  
  first_column <- colnames(input)[1]
  
  iterations <- unique(input[[first_column]])
  
  output <- tibble(!!sym(first_column) := iterations, func_act = NA, p_value = NA, p_adj = NA, 
                   rel_abd_prot = NA, rel_abd_gene = NA)
  
  for (i in 1:length(iterations)) {
    iteration <- iterations[i]
    
    comparison <- input %>%
      filter(!!sym(first_column) == iteration)
    
    wilcox <- wilcox.test(comparison$rel_abd_gene, comparison$rel_abd_prot, paired = T)
    
    p <- wilcox$p.value
    
    output$p_value[i] <- p
    
    output$func_act[i] <- median(comparison$func_act) # changed to median
    output$rel_abd_prot[i] <- median(comparison$rel_abd_prot)
    output$rel_abd_gene[i] <- median(comparison$rel_abd_gene)
    
  }
  
  # adjust p_value
  
  output$p_adj <- p.adjust(output$p_value, method = "holm")
  
  output <- output %>%
    mutate(sig = ifelse(p_adj < 0.05, "sig", "not sig")) %>%
    mutate(rank_asc = row_number(func_act),
           rank_desc = row_number(-func_act),
           label = ifelse(sig == "sig" & (rank_desc <= 10 | rank_asc <= 10), !!sym(first_column), ""))
  
  p <- ggplot(output, aes(x = func_act, y = -log10(p_value))) +
    geom_vline(xintercept = 0, color = "grey20") +
    geom_hline(yintercept = 0, color = "grey20") + 
    geom_point(aes(size = rel_abd_prot), show.legend = T, color = "#023743FF", alpha = .7) +
    geom_point(aes(size = rel_abd_gene), show.legend = F, color = "#FED789FF") +
    geom_text_repel(aes(label = label), color = "black", show.legend = F, max.overlaps = Inf, box.padding = 0.5,
                    min.segment.length = 0, max.time = 3, size = 3, nudge_y = -.1) +
    scale_color_manual(values = c("grey20", "red")) +
    #scale_size_manual(values = c(1,2)) +
    labs(x = "log2 fold change", y = "-log10 P-value", size = "rel. abd. %") 
  
  print(p)
  
  return(output)
}

# functional distance between metagenomic and metaproteomic frame
# func dist function

calculate_func_dist <- function(integrated_df) {
  
  func_dist <- tibble(sampleid = unique(integrated_df$sampleid),
                      func_dist = NA)
  
  for (i in 1:nrow(func_dist)) {
    current_sampleid <- func_dist$sampleid[i]
    
    filtered_matrix <- integrated_df %>%
      filter(sampleid == current_sampleid) %>%
      dplyr::select(rel_abd_gene, rel_abd_prot) %>%
      as.matrix() %>%
      t()
    
    bray <- vegdist(filtered_matrix)
    
    func_dist$func_dist[i] <- bray
    
  }
  
  func_dist_out <- func_dist %>%
    inner_join(meta, by = "sampleid")
  
  p <- ggplot(func_dist_out, aes(x = diet, y = func_dist))+
    geom_boxplot(outliers = F) +
    geom_quasirandom(aes(color = animal)) +
    facet_grid(~matrix) +
    scale_color_manual(values = colors)
  print(p)
  
  return(func_dist_out)
}

########################################
# Functions to save top features
########################################



save_top_features <- function(df, selected_matrix, top_features = 10, save, save_name) {
  feature_column <- colnames(df)[1]
  abundance_column <- colnames(df)[3]
  
  print(paste("Feature column is", feature_column))
  print(paste("Abundance column is called", abundance_column))
  
  if (selected_matrix == "ileal digesta") {
    filtered_df <- filter_ileum(df)
  } else if (selected_matrix == "faeces") {
    filtered_df <- filter_faeces(df)
  }
  
  output_df <- filtered_df %>%
    dplyr::select(!!sym(feature_column), sampleid, !!sym(abundance_column)) %>%
    group_by(!!sym(feature_column)) %>%
    mutate(mean_rel_abd = mean(!!sym(abundance_column))) %>%
    ungroup() %>%
    mutate(ranking = dense_rank(desc(mean_rel_abd))) %>%
    filter(ranking <= top_features) %>%
    dplyr::select(-c(mean_rel_abd, ranking))
  # save if wanted
  if (isTRUE(save)) {
    matrix_abbr <- ifelse(selected_matrix == "ileal digesta", "il", "fa")
    top_abbr <- ifelse(top_features < Inf, "_top_", "_all_")
    saveRDS(output_df, file = paste0("clean/", get_script_number(), top_abbr, matrix_abbr, "_", save_name, ".RDS"))
  }
  return(output_df)
}

save_top_features_k_table <- function(df, selected_matrix, top_features = 10, save, save_name) {
  feature_column <- colnames(df)[1]
  abundance_column <- colnames(df)[4]
  
  print(paste("Feature column is", feature_column))
  print(paste("Abundance column is called", abundance_column))
  
  if (selected_matrix == "ileal digesta") {
    filtered_df <- filter_ileum(df)
  } else if (selected_matrix == "faeces") {
    filtered_df <- filter_faeces(df)
  }
  
  output_df <- filtered_df %>%
    dplyr::select(!!sym(feature_column), rank, sampleid, !!sym(abundance_column)) %>%
    filter(rank %in% c("F", "G", "S")) %>%
    group_by(!!sym(feature_column), rank) %>%
    mutate(mean_rel_abd = mean(!!sym(abundance_column))) %>%
    ungroup() %>%
    group_by(rank) %>% # ranking within ranks
    mutate(ranking = dense_rank(desc(mean_rel_abd))) %>%
    ungroup() %>%
    filter(ranking <= top_features) %>%
    dplyr::select(-c(mean_rel_abd, ranking, rank))
  # save if wanted
  if (isTRUE(save)) {
    matrix_abbr <- ifelse(selected_matrix == "ileal digesta", "il", "fa")
    top_abbr <- ifelse(top_features < Inf, "_top_", "_all_")
    saveRDS(output_df, file = paste0("clean/", get_script_number(), top_abbr, matrix_abbr, "_", save_name, ".RDS"))
  }
  return(output_df)
}


##################
# Functional reduncancy
####################

# Function that calculates functional redundancy
# takes a table with proteinid, respective taxonomic classification (G), functional classification (kegg_cog), sampleid, and relative abundance
# takes table (on genus level) with name, sampleid, and relative abundance
# meta with sampleid needed

calculate_functional_redundancy <- function(pro_gen_cog, genus_table, meta) {
  # build PCN
  
  for (i in unique(meta$sampleid)) {
    #print(i)
    pro_gen_cog_sam_sum_wide <- pro_gen_cog %>%
      dplyr::select(G, kegg_cog, sampleid, rel_abd) %>%
      filter(sampleid == i) %>%
      filter(!is.na(kegg_cog)) %>% 
      group_by(G, kegg_cog) %>%
      summarise(sum = sum(rel_abd)) %>%
      pivot_wider(names_from = G, values_from = sum)
    assign(paste0("pcn", i), pro_gen_cog_sam_sum_wide)
  }
  
  # build 01 PCNs
  # 
  # for (i in unique(meta$sampleid)) {
  #   print(i)
  #   PCN_table_single_filtered_sum_col <- get(paste0("pcn", i)) %>%
  #     replace(is.na(.), 0) %>%
  #     pivot_longer(-kegg_cog, names_to = "G", values_to = "value") %>%
  #     group_by(G) %>%
  #     filter(sum(value) > 0) %>%
  #     ungroup() %>%
  #     group_by(kegg_cog) %>%
  #     filter(sum(value) > 0) %>%
  #     ungroup() %>%
  #     mutate(value = ifelse(value == 0, 0, 1)) %>%
  #     pivot_wider(names_from = kegg_cog, values_from = value)
  #   assign(paste0("pcn01", i), PCN_table_single_filtered_sum_col)
  # }
  
  # dij calculation
  
  dij <- function(mat, x, y) {
    min <- sum(pmin(mat[which(rownames(mat) == x),], mat[which(rownames(mat) == y),]))
    max <- sum(pmax(mat[which(rownames(mat) == x),], mat[which(rownames(mat) == y),]))
    dij <- 1 - (min / max)
    return(dij)
  }
  
  
  pb <- txtProgressBar(min = 1, max = length(unique(meta$sampleid)), style = 3) # initialize progessbar
  for (idx in 1:length(unique(meta$sampleid))) {
    i <- unique(meta$sampleid)[idx]
    setTxtProgressBar(pb, idx) # print progessbar
    PCN <- get(paste0("pcn", i)) %>%
      replace(is.na(.), 0) %>%    
      pivot_longer(-kegg_cog, names_to = "G", values_to = "value") %>%
      group_by(kegg_cog) %>%
      filter(sum(value) > 0) %>%
      ungroup() %>%
      group_by(G) %>%
      filter(sum(value) > 0) %>%
      mutate(norm = value / sum(value)) %>%
      ungroup() %>%
      dplyr::select(-value) %>%
      pivot_wider(names_from = kegg_cog, values_from = norm)
    PCNmat <- as.matrix(PCN[,-1])
    rownames(PCNmat) <- PCN$G
    dijmat <- matrix(nrow = nrow(PCNmat), ncol = nrow(PCNmat))
    rownames(dijmat) <- PCN$G
    colnames(dijmat) <- PCN$G
    
    for (j in rownames(PCNmat)) {
      #print(j)
      for (k in rownames(PCNmat)) {
        #print(k)
        dijmat[j,k] <- dij(mat = PCNmat, x = j, y = k)
      }
    }
    assign(paste0("dij", i), dijmat)
  }
  close(pb) # close progessbar
  
  # FR calculation
  
  genus_table_p_matrix <- genus_table %>%
    group_by(sampleid) %>%
    mutate(rel_abd = rel_abd/sum(rel_abd)) %>%
    pivot_wider(names_from = sampleid, values_from = rel_abd)
  
  output <- matrix(ncol = 4, nrow = length(unique(meta$sampleid))) # matrix to store output
  rownames(output) <- unique(meta$sampleid)
  colnames(output) <- c("FR", "nFR", "GSI", "FD")
  
  for (i in unique(meta$sampleid)) {
    dijmat <- get(paste0("dij", i))
    genus_table <- genus_table_p_matrix %>% # discarding genera not in dij
      filter(name %in% rownames(dijmat)) %>%
      mutate(name = factor(name, levels = rownames(dijmat))) %>%
      arrange(name)
    for (z in 1:length(genus_table$name)) {
      if (genus_table$name[z] != rownames(dijmat)[z]) {
        print("Names do not match!")
      }
    }
    genus_table$pipj <- 0
    genus_table$dijpipj <- 0
    
    for (j in row_number(genus_table)) {
      genus_table$pipj[j] <- 
        sum(genus_table[[i]][j] * 
              genus_table[[i]][-j])
      
      genus_table$dijpipj[j] <- 
        sum(genus_table[[i]][j] * genus_table[[i]][-j] * 
              (1 - dijmat[which(rownames(dijmat) == genus_table$name[j]),
                          which(rownames(dijmat) != genus_table$name[j])]))
    }
    GSI <- sum(genus_table$pipj) 
    FR <- sum(genus_table$dijpipj)
    output[which(rownames(output) == i), "FR"] <- FR
    output[which(rownames(output) == i), "nFR"] <- FR / GSI
    output[which(rownames(output) == i), "GSI"] <- GSI
    output[which(rownames(output) == i), "FD"] <- GSI - FR
  }
  
  return(output)
}
