library(SmCCNet)
library(mixOmics)
library(RCy3)
library(igraph)
library(KEGGREST)
library(pheatmap)

source(here("30_omics_functions.R"))

# filter function for binary outcome

filter_input_smccn <- function(omics_list, meta, include) {
  meta_filtered <- meta %>%
    filter(diet %in% include) %>%
    mutate(diet = ifelse(diet == include[1], "0", "1")) %>%
    mutate(diet = as.factor(diet)) #as factor for PLS
  
  list_filtered <- omics_list
  for (i in 1:length(list_filtered)) {
    list_filtered[[i]] <- list_filtered[[i]][rownames(list_filtered[[i]]) %in% meta_filtered$sampleno,]
  }
  
  output_list <- list()
  output_list[[1]] <- list_filtered
  output_list[[2]] <- meta_filtered
  return(output_list)
}

# function to load a subnetwork from the current save_folder

load_subnetwork <- function(network) {
  subnetworks <- list.files(save_folder)
  network_file <- subnetworks[which(endsWith(subnetworks, str_c("_", network, ".Rdata")))]
  print(paste("load", network_file))
  load(file.path(save_folder, network_file), envir = globalenv())
}

# load the all keggs from the correlation matrix and remove the prefix

load_keggs <- function() {
  features <- rownames(correlation_sub)
  keggs <- features %>%
    str_remove("[a-z]+_") 
  keggs <- keggs[which(str_detect(keggs, "K\\d{5}"))]
  cat(keggs)
  return(keggs)
}

# do enrichment analysis with the kegg_kos and the universe

enrich_keggs <- function(kegg_kos) {
  enrichment <- enrichKO(kegg_kos, universe = universe)
  if (length(which(enrichment@result$p.adjust<enrichment@pvalueCutoff)) > 0) {
    p <- enrichplot::dotplot(enrichment)
    print(p)
  } else {
    print(paste("No enriched KEGGs found"))
  }
  return(enrichment)
}

# add information to kegg_kos

annotate_keggs <- function(kegg_kos) {
  out_list <- list()
  for (i in 1:length(kegg_kos)) {
    ko <- kegg_kos[i]
    print(ko)
    
    out <- tryCatch(
      {
        kegg_info <- keggGet(ko)[[1]]
        tibble(kegg_ko = ko,
               name = kegg_info$NAME,
               pathways = str_c(kegg_info$PATHWAY, collapse = ","),
               pathway_ids = str_c(names(kegg_info$PATHWAY), collapse = ","),
               modules = str_c(kegg_info$MODULE, collapse = ","),
               module_ids = str_c(names(kegg_info$MODULE),collapse= ","),
               reactions = str_c(kegg_info$REACTION, collapse = ","),
               reaction_ids = str_c(names(kegg_info$REACTION), collapse = ","))
      },
      error = function(e) {
        tibble(kegg_ko = ko,
               name = NA,
               pathways = NA,
               pathway_ids = NA,
               modules = NA,
               module_ids = NA,
               reactions = NA,
               reaction_ids = NA)
      }
    )
    
    out_list[[i]] <- out
  }
  out <- do.call("rbind", out_list)
  
  return(out)
}

#############
# diablo functions
##################

# detect pairwise correlations
diablo_pairwise_correlations <- function(omics_list) {
  vector <- c()
  for (i in 1:length(omics_list)) {
    X <- omics_list[[i]]
    for (j in 1:length(omics_list)) {
      Y <- omics_list[[j]]
      if (j > i) {
        print(paste(names(omics_list)[i], "vs", names(omics_list)[j]))
        pls <- spls(X, Y, keepX = c(min(ncol(X), 20),min(ncol(X), 20)), 
                    keepY = c(min(ncol(Y), 20),min(ncol(Y), 20)))
        plotVar(pls, var.names = T, cutoff = 0.5)
        print(cor(pls$variates$X, pls$variates$Y))
        vector[i] <- cor(pls$variates$X, pls$variates$Y)[1,1]
        Sys.sleep(2)
      }
    }
  }
  
  print(paste("Mean correlation:", mean(vector)))
  return(mean(vector))
}

# set up and test model, return final model

diablo_auto_model <- function(omics_list, meta, design_value, ncomp_value = NULL) {
  # measure execution time
  start_time <- Sys.time()
  
  #create design matrix
  design = matrix(design_value, ncol = length(omics_list), nrow = length(omics_list), 
                  dimnames = list(names(omics_list), names(omics_list)))
  diag(design) = 0 # set diagonal to 0s
  
  # initialise 
  diablo_init <- block.splsda(X = omics_list, Y = meta, ncomp = 10, design = design)
  
  # test for best number of comps
  diablo_perf <- perf(diablo_init, validation = "loo",
                      progressBar = T, seed = 1112, nrepeat = 5) # nrepeat to obtain optimal ncomp numbers
  par(mfrow = c(1,1))
  plot(diablo_perf)
  
  print(diablo_perf$choice.ncomp)
  
  if (is.null(ncomp_value)) {
    #ncomp = min(unlist(diablo_perf$choice.ncomp))
    ncomp = as.integer(readline(prompt = "Enter number of components"))
    print(paste("Choice of components:", ncomp))
  } else {
    ncomp = ncomp_value
  }

  
  # test different numbers of features
  test.keepX <- list()
  for (i in 1:length(omics_list)) {
    test.keepX[[names(omics_list)[i]]] = seq(5, min(ncol(omics_list[[i]]), 45), 10)
  }
  
  diablo_tune <- tune.block.splsda(X = omics_list, Y = meta,
                                   ncomp = ncomp, test.keepX = test.keepX,
                                   design = design, validation = "loo",
                                   dist = "mahalanobis.dist",
                                   BPPARAM = BiocParallel::SnowParam(workers = 8), # 8 threads
                                   progressBar = T, seed = 1112)
  plot(diablo_tune)
  
  list.keepX <- diablo_tune$choice.keepX
  print("Choice of features:")
  print(list.keepX)
  
  # final model
  
  diablo_final <- block.splsda(X = omics_list, Y = meta, ncomp = ncomp,
                               keepX = list.keepX, design = design)
  
  # measure execution time
  end_time <- Sys.time()
  print(end_time-start_time)
  
  return(diablo_final)
}

# Create every plot with final diablo model

plot_diablo <- function(final_model, interactive_network = FALSE, cutoff_value = 0.8, save_name = NULL, 
                        return_network = FALSE, return_full_matrix = FALSE) {
  plotDiablo(final_model, ncomp = min(final_model$ncomp[1], 2))
  if (final_model$ncomp[1] > 1) {
    plotIndiv(final_model)
    plotArrow(final_model)
    plotVar(final_model)
  }

  par(mfrow = c(1,1))
  
  # helper to select use components for network
  comps <- names(final_model$ncomp)[1:(length(final_model$ncomp)-1)]
  list_comps <- list()
  for (i in 1:length(comps)) {
    list_comps[[i]] <- seq(1, final_model$ncomp[1])
    names(list_comps)[i] <- comps[i]
  }
  
  matrix <- circosPlot(final_model, cutoff = cutoff_value, size.variables = 1, line = T, size.labels = 1.5)
  network <- network(final_model, cutoff = cutoff_value, blocks = seq(1, length(final_model$names$blocks)-1),
              comp = list_comps, interactive = interactive_network,
          size.node = 0.1, cex.node.name = 0.5, lwd.edge = 1, symkey = T,
          color.edge = c(color.jet(100)[1:25], color.jet(100)[76:100]),
          save = if (is.null(save_name)) save_name else "jpeg",
          name.save = if (is.null(save_name)) NULL else { 
                             paste0("plots/", get_script_number(), "_diablo_network_", save_name)})
  
  if (!isTRUE(return_full_matrix)) {
    diag(matrix) <- diag(matrix) / 100 # helper to mitigate filtering due to self correlation
    matrix[which(abs(matrix)<cutoff_value)] <- 0
    matrix <- matrix[which(rowSums(matrix) != 0),
                     which(colSums(matrix) != 0)]
    diag(matrix) <- diag(matrix) * 100
  }

  pheatmap(matrix, col = color.jet(100), scale = "none", fontsize = 8)

  plotLoadings(final_model, comp = 1)
  plotLoadings(final_model, comp = 2)
  cimDiablo(final_model, legend.position = "left", transpose = T)
  if (isTRUE(return_network)) {
    return(network)
  } else if (isTRUE(return_full_matrix)) {
    return(circosPlot(final_model, cutoff = 0))
  } else {
    return(matrix)
  }

}

# filter correlation matrix to a subnetwork given a vector of names and create a heatmap

filter_submatrix <- function(matrix, vector, save_name = NULL) {
  not_found <- vector[!vector %in% rownames(matrix)]
  if (length(not_found) > 0){
    print(paste("Not found:", not_found))
  }

  index_vector <- which(rownames(matrix) %in% vector)
  submatrix <- matrix[index_vector, index_vector]
  if (!is.null(save_name)) {
    filename = paste0("plots/", get_script_number(), "_heatmap_", save_name, ".jpeg")
  } else {
    filename = NA
  }
  pheatmap(submatrix, filename = filename,
           angle_col = 45, width = 10)
  return(submatrix)
}

# create a igraph from a adjacency matrix

create_igraph_from_matrix <- function(matrix, save_name = NULL) {
  g <- graph_from_adjacency_matrix(matrix, mode = "undirected", weighted = T, diag = F)
  layout_fr <- layout_with_fr(g, weights = 1-(abs(E(g)$weight)))
  # assign colors to weight values
  color_index <- round((E(g)$weight + 1) / 2 * 99) + 1
  E(g)$color <- color.jet(100)[color_index]
  E(g)$width <- .5
  node_color <- c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")
  color_vector <- case_when(str_detect(V(g)$name, "^il|^fa|^pc|^hg|^enz") ~ 1,
                            str_detect(V(g)$name, "^g_") ~ 2,
                            str_detect(V(g)$name, "^p_") ~ 3,
                            str_detect(V(g)$name, "^ssc_") ~ 4,
                            str_detect(V(g)$name, "_PEA$") ~ 5,
                            str_detect(V(g)$name, "^m_") ~ 6,
                            .default = 0)
  V(g)$color <- node_color[color_vector]
  V(g)$label.family <- "sans"
  V(g)$label.font <- 1
  V(g)$label.color <- "black"
  V(g)$label.cex <- .15
  if (!is.null(save_name)) jpeg(filename = paste0("plots/", get_script_number(), "_igraph_", save_name, ".jpeg"),
                                width = 4000, height = 4000, unit="px", res = 1200, pointsize = 20)
  par(mar = c(0,0,0,0))
  plot(g, layout = layout_fr, margin = c(0,0,0,0), rescale = TRUE,
       vertex.shape = "rectangle", 
       vertex.size = str_width(V(g)$name)*2,
       vertex.size2 = 4,
       vertex.frame.width = .1)
  legend("topright",
         legend = c(1,NA,NA,NA,NA, 0, NA,NA,NA,NA,-1),
         fill = rev(color.jet(11)),
         border = NA,
         y.intersp = 0.5,
         cex = .2, text.font = 1)
  legend("bottomright",
         legend = c("Nutrition", "Metagenomics", "Metaproteomics", "Host proteins", "Pea proteins", "Metabolomics"),
         fill = node_color,
         cex = .2)
  if (!is.null(save_name)) dev.off() 
}

# takes a list with feature_x, feature_y and correlation_type column and creates a adjacency matrix

list_to_matrix <- function(list) {
  matrix1 <- list %>%
    mutate(value = ifelse(correlation_type == "+", 1, -1)) %>%
    dplyr::select(x = feature_x, y = feature_y, value)
  matrix2 <- rbind(dplyr::select(matrix1, feature_x = x, feature_y = y, value), 
                   dplyr::select(matrix1, feature_x = y, feature_y = x, value)) %>%
    pivot_wider(names_from = feature_y, values_from = value) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.x), 0, .x))) %>%
    column_to_rownames("feature_x") %>%
    as.matrix()
  matrix2 <- matrix2[, sort(colnames(matrix2))]
  matrix2 <- matrix2[sort(rownames(matrix2)),]
  return(matrix2)
}


