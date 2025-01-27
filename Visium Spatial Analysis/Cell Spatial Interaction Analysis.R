## Previously developed in Mascharak & Guo et al. Cell Reports Medicine 2023.

## Initialization
library(AdjacencyScore)
library(igraph)
library(ggraph)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
source("Additional Functions.R")


# Prior to this point, analyze Visium data and integrate with scRNA-seq reference
# using the standard Seurat pipeline (stored in object "vsm").
# See official vignette at https://satijalab.org/seurat/articles/spatial_vignette.html.


## Pre-Analysis
# Load in sample-level metadata (.csv should have one row per sample, one column titled "Sample.Name" for the aligned sample ID, 
# and one column titled "Group" with the group assignment of each sample. Groups should be titled as GroupName1.Group, GroupName2.Group, etc.
# A header column, representing the slice IDs, and a header row should also be included.)
loaded_metadata <- read.delim(file = "sample_metadata.csv", header = TRUE, row.names=1, sep = ",")
loaded_metadata <- t(loaded_metadata)
# Store experimental group in Seurat
temp_metadata <- loaded_metadata
colnames(temp_metadata) <- loaded_metadata["Sample.Name",]
vsm$group <- temp_metadata["Group",vsm$sample.name]
rm(temp_metadata)
# Set default assay to cell type probabilities
DefaultAssay(vsm) <- "predictions"
# Convert image coordinates to numeric format 
for (sample in names(vsm@images))
{
  for (coordinate in names(vsm@images[[sample]]@coordinates))
  {
    vsm@images[[sample]]@coordinates[[coordinate]] <- as.numeric(vsm@images[[sample]]@coordinates[[coordinate]])
  }
}
# Store x/y coordinates of spots
cell_locations <- list()
for (sample in names(vsm@images))
{
  cell_locations[[sample]] <- as.matrix(vsm@images[[sample]]@coordinates[,c("row","col")])
}
# Store cell probability matrices
cell_probability_matrix <- list()
for (sample in names(vsm@images))
{
  cell_probability_matrix[[sample]] <- t(as.matrix(vsm@assays$predictions@data[rownames(vsm@assays$predictions@data) != "max", vsm$sample.name == loaded_metadata["Sample.Name",sample]]))
}
# Specify discrete categories of interest
discrete_comparisons <- c("Group")
discrete_groups <- unique(loaded_metadata[discrete_comparisons,])
discrete_groups <- discrete_groups[!is.na(discrete_groups)]
# Specify pairwise comparisons for discrete categories (right side will be colored red in later plots)
discrete_pairs <- c("D0.Group vs. D14.Group","D0.Group vs. D35.Group","D14.Group vs. D35.Group")


## Spatial Analysis - Sample Level
# Specify whether to output individual vs. sample-level scores, normalize
# interaction scores by number of spots, and/or log-scale scores
individual_setting <- FALSE
normalize_setting <- TRUE
log_setting <- TRUE
# Generate all possible cell pairs
cell_types <- sort(colnames(cell_probability_matrix[[1]]))
cell_types <- cell_types[cell_types != "Red Blood Cells"] # Exclude red blood cells from downstream analysis
# cell_types <- cell_types[grepl("ATI|Fibroblast",cell_types)] # Calculate interactions for fibroblasts and ATI/ATII cells only
cell_pairs <- t(combn(cell_types,2))
for (cell in cell_types) 
{
  cell_pairs <- rbind(cell_pairs, c(cell,cell))
}
# Create empty matrix to store for excessively small samples 
empty_matrix <- matrix(nrow = nrow(cell_pairs), ncol = 3)
rownames(empty_matrix) <- paste(cell_pairs[,1], cell_pairs[,2], sep = " - ")
colnames(empty_matrix) <- c("f1","f2","score") 
empty_matrix[,"f1"] <- cell_pairs[,1]
empty_matrix[,"f2"] <- cell_pairs[,2]
empty_matrix[rownames(empty_matrix),"score"] <- 0 
# Generate interaction scores for a variety of k parameter values
# k=6 was utilized in the present study to assess adjacent spots (i,e, local neighborhoods) only
k_interaction_scores <- list()
k_interaction_network <- list()
for (interaction_k in c(6,18,36,60,90))
{
  print(paste0("BEGINNING CALCULATIONS FOR k=",interaction_k,"."))
  # Generate interaction scores and network for each sample 
  interaction_scores <- list()
  interaction_network <- list()
  for (sample in names(cell_locations))
  {
    print(paste0("Calculating interaction scores for sample: ",sample,"."))
    # If too few cells, populate with empty matrix
    if (nrow(cell_locations[[sample]]) < interaction_k)
    {
      interaction_scores[[sample]] <- empty_matrix
    }
    else
    {
      # Generate adjacency matrix based on cell locations
      adjacency_matrix <- knn_graph(cell_locations[[sample]], k = interaction_k)
      # Calculate interaction scores as dot-product of paired cell feature vectors with adjacency matrix (f1*j*f2)
      interaction_scores[[sample]] <- CalculateInteractionScores(adjacency_matrix = adjacency_matrix, cell_matrix = t(cell_probability_matrix[[sample]]), cell_pairs = cell_pairs, calculate_individual_scores = individual_setting, normalize = normalize_setting, log_scale = log_setting)
      # Generate interaction matrix and resulting network
      interaction_score_matrix <- AdjScoreHeatmapMatrix(interaction_scores[[sample]])
      class(interaction_score_matrix) <- "numeric"
      interaction_network[[sample]] <- graph_from_adjacency_matrix(interaction_score_matrix, weighted = TRUE)
    }
  }
  # Store interaction scores and networks
  k_interaction_scores[[as.character(interaction_k)]] <- interaction_scores
  k_interaction_network[[as.character(interaction_k)]] <- interaction_network
  rm(adjacency_matrix,interaction_score_matrix)
}
rm(interaction_scores,interaction_network)


## Calculation of Network Centrality
# Load networks
interaction_k <- 6
interaction_network <- k_interaction_network[[as.character(interaction_k)]]
# Calculate weighted eigenvector centrality
ev_centrality <- matrix(nrow = length(cell_types), ncol = 6)
rownames(ev_centrality) <- cell_types
colnames(ev_centrality) <- c("D0.Group","D14.Group","D35.Group","D0.Group.SD","D14.Group.SD","D35.Group.SD")
for (group in colnames(ev_centrality)[1:3])
{
  relevant_samples <- colnames(loaded_metadata[,loaded_metadata["Group",] == group])
  relevant_scores <- matrix(nrow = nrow(ev_centrality), ncol = length(relevant_samples))
  rownames(relevant_scores) <- rownames(ev_centrality)
  colnames(relevant_scores) <- relevant_samples
  for (sample in relevant_samples)
  {
    relevant_scores[,sample] <- as.numeric(evcent(interaction_network[[sample]])$vector[rownames(relevant_scores)])
  }
  for (node in rownames(ev_centrality))
  {
    ev_centrality[node,group] <- mean(relevant_scores[node,], na.rm = TRUE)
    ev_centrality[node,paste0(group,".SD")] <- sd(relevant_scores[node,], na.rm = TRUE)
  }
}


## Spatial Analysis - Group Level
# Calculate average interaction scores for each group and produce group-level interaction network
k_group_interaction_scores <- list()
k_group_interaction_scores_with_zeroes <- list()
k_group_interaction_network <- list()
for (interaction_k in c(6,18,36,60,90))
{
  # Load variables for current k parameter
  interaction_scores <- k_interaction_scores[[as.character(interaction_k)]]
  # Calculate average interaction scores for each group and produce resulting interaction network
  group_interaction_scores <- list()
  group_interaction_scores_with_zeroes <- list()
  group_interaction_network <- list()
  for (group in discrete_groups)
  {
    # Define group-specific samples and variable for storing group-specific interaction scores
    metadata_category <- strsplit(group, split = "\\.")[[1]][2]
    group_samples <- names(which(loaded_metadata[metadata_category,] == group))
    group_samples <- group_samples[group_samples %in% names(cell_locations)]
    group_interaction_scores[[group]] <- matrix(nrow = nrow(cell_pairs), ncol=3)
    colnames(group_interaction_scores[[group]]) <- c("f1","f2","score") 
    group_interaction_scores[[group]][,"f1"] <- cell_pairs[,1]
    group_interaction_scores[[group]][,"f2"] <- cell_pairs[,2]
    rownames(group_interaction_scores[[group]]) <- paste(group_interaction_scores[[group]][,"f1"], group_interaction_scores[[group]][,"f2"], sep = " - ")
    # For each possible cell pair, calculate the average interaction score across all group-specific samples
    for (pair in 1:nrow(cell_pairs))
    {
      score_sum <- 0
      for (sample in group_samples)
      {
        matching_row_1 <- (interaction_scores[[sample]][,"f1"] == cell_pairs[pair,1]) & (interaction_scores[[sample]][,"f2"] == cell_pairs[pair,2])
        matching_row_2 <- (interaction_scores[[sample]][,"f2"] == cell_pairs[pair,1]) & (interaction_scores[[sample]][,"f1"] == cell_pairs[pair,2])
        if (sum(matching_row_1) == 1)
        {
          if (!is.na(interaction_scores[[sample]][matching_row_1,"score"]))
          {
            score_sum <- score_sum + as.numeric(interaction_scores[[sample]][matching_row_1,"score"])
          }
        }
        else if (sum(matching_row_2) == 1)
        {
          if (!is.na(interaction_scores[[sample]][matching_row_2,"score"]))
          {
            score_sum <- score_sum + as.numeric(interaction_scores[[sample]][matching_row_2,"score"])
          }
        }
      }
      if (score_sum != 0)
      {
        group_interaction_scores[[group]][pair,"score"] <- score_sum/length(group_samples)
      }
    }
    # Save full interaction score variable (including zero/NA values) for later usage in differential analysis
    group_interaction_scores_with_zeroes[[group]] <- group_interaction_scores[[group]]
    group_interaction_scores_with_zeroes[[group]][is.na(group_interaction_scores_with_zeroes[[group]][,"score"]),"score"] <- 0
    # Remove non-existing interactions from group-averaged interaction score variable
    group_interaction_scores[[group]] <- group_interaction_scores[[group]][!is.na(group_interaction_scores[[group]][,"score"]),]
    # Generate interaction matrix and resulting network for plotting
    group_interaction_score_matrix <- AdjScoreHeatmapMatrix(group_interaction_scores[[group]])
    class(group_interaction_score_matrix) <- "numeric"
    group_interaction_network[[group]] <- graph_from_adjacency_matrix(group_interaction_score_matrix, weighted = TRUE)
  }
  # Store group interaction scores and networks
  k_group_interaction_scores[[as.character(interaction_k)]] <- group_interaction_scores
  k_group_interaction_scores_with_zeroes[[as.character(interaction_k)]] <- group_interaction_scores_with_zeroes
  k_group_interaction_network[[as.character(interaction_k)]] <- group_interaction_network
  rm(group_interaction_score_matrix)
}
rm(group_interaction_scores,group_interaction_scores_with_zeroes,group_interaction_network)


## Spatial Analysis - Differential Level
# Calculate differential interaction scores and produce resulting differential interaction network
k_differential_interaction_scores <- list()
k_differential_interaction_network <- list()
k_differential_interaction_color <- list()
for (interaction_k in c(6,18,36,60,90))
{
  # Load variables for current k parameter
  interaction_scores <- k_interaction_scores[[as.character(interaction_k)]]
  group_interaction_scores_with_zeroes <- k_group_interaction_scores_with_zeroes[[as.character(interaction_k)]]
  # For discrete groups, calculate differential interaction scores and generate resulting differential interaction networks
  differential_interaction_scores <- list()
  differential_interaction_network <- list()
  differential_interaction_color <- list()
  for (comparison in discrete_pairs)
  {
    # Calculate differences in interaction scores between groups
    group_1 <- str_split(comparison," vs. ")[[1]][1]
    group_2 <- str_split(comparison," vs. ")[[1]][2]
    differential_interaction_scores[[comparison]] <- group_interaction_scores_with_zeroes[[group_1]]
    differential_interaction_scores[[comparison]][,"score"] <- as.numeric(group_interaction_scores_with_zeroes[[group_2]][,"score"]) - as.numeric(group_interaction_scores_with_zeroes[[group_1]][,"score"])
    # Remove non-existing interactions from group-averaged interaction score variable
    differential_interaction_scores[[comparison]] <- differential_interaction_scores[[comparison]][differential_interaction_scores[[comparison]][,"score"] != 0,]
    # Generate differential interaction score matrix
    differential_interaction_score_matrix <- AdjScoreHeatmapMatrix(differential_interaction_scores[[comparison]])
    class(differential_interaction_score_matrix) <- "numeric"
    # Convert differential interaction score matrix to absolute values and separate color variables for plotting
    differential_interaction_color[[comparison]] <- differential_interaction_score_matrix
    differential_interaction_color[[comparison]][differential_interaction_score_matrix < 0] <- "Negative"
    differential_interaction_color[[comparison]][differential_interaction_score_matrix > 0] <- "Positive"
    differential_interaction_score_matrix <- abs(differential_interaction_score_matrix)
    # Generate differential interaction network
    differential_interaction_network[[comparison]] <- graph_from_adjacency_matrix(differential_interaction_score_matrix, weighted = TRUE)
  }
  # Store differential interaction scores, networks, and colors
  k_differential_interaction_scores[[as.character(interaction_k)]] <- differential_interaction_scores
  k_differential_interaction_network[[as.character(interaction_k)]] <- differential_interaction_network
  k_differential_interaction_color[[as.character(interaction_k)]] <- differential_interaction_color
  rm(differential_interaction_score_matrix)
}
rm(differential_interaction_scores,differential_interaction_network,differential_interaction_color)


## Spatial Analysis - Individual Interaction Level
# Calculate group-level summaries of individual interaction scores
k_interaction_score_overview <- list()
k_group_interaction_score_overview <- list()
k_group_interaction_score_overview_pvalues <- list()
for (interaction_k in c(6,18,36,60,90))
{  
  # Load variables for current k parameter
  interaction_scores <- k_interaction_scores[[as.character(interaction_k)]]
  # Extract individual interaction scores per sample
  interaction_score_overview <- matrix(nrow = nrow(cell_pairs), ncol = length(names(cell_locations)))
  rownames(interaction_score_overview) <- paste(cell_pairs[,1],cell_pairs[,2], sep = " - ")
  colnames(interaction_score_overview) <- names(cell_locations)
  for (interaction in rownames(interaction_score_overview))
  {
    for (sample in colnames(interaction_score_overview))
    {
      interaction_score_overview[interaction,sample] <- as.numeric(interaction_scores[[sample]][interaction,"score"])
    }
  }
  # Calculate mean and sd of individual interactions per group
  group_interaction_score_overview <- matrix(nrow = length(discrete_groups), ncol = nrow(interaction_score_overview)*2+1)
  rownames(group_interaction_score_overview) <- discrete_groups
  colnames(group_interaction_score_overview) <- c("Group",rownames(interaction_score_overview),paste0(rownames(interaction_score_overview),".sd"))
  group_interaction_score_overview[,"Group"] <- rownames(group_interaction_score_overview)
  for (group in rownames(group_interaction_score_overview))
  {
    metadata_category <- strsplit(group, split = "\\.")[[1]][2]
    group_samples <- names(which(loaded_metadata[metadata_category,] == group))
    group_samples <- group_samples[group_samples %in% names(cell_locations)]
    for (interaction in rownames(interaction_score_overview))
    {
      group_interaction_score_overview[group,interaction] <- mean(interaction_score_overview[interaction,group_samples])
      group_interaction_score_overview[group,paste0(interaction,".sd")] <- sd(interaction_score_overview[interaction,group_samples])
    }
  }
  group_interaction_score_overview <- as.data.frame(group_interaction_score_overview)
  group_interaction_score_overview$Group <- factor(group_interaction_score_overview$Group, levels = group_interaction_score_overview$Group)
  # Calculate p values of group differences in interaction score
  group_interaction_score_overview_pvalues <- matrix(nrow = nrow(interaction_score_overview), ncol = length(discrete_pairs))
  rownames(group_interaction_score_overview_pvalues) <- rownames(interaction_score_overview)
  colnames(group_interaction_score_overview_pvalues) <- discrete_pairs
  for (interaction in rownames(interaction_correlations))
  {
    for (comparison in discrete_pairs)
    {
      group_1 <- str_split(comparison," vs. ")[[1]][1]
      group_2 <- str_split(comparison," vs. ")[[1]][2]
      metadata_category <- strsplit(group_1, split = "\\.")[[1]][2]
      group_1_samples <- names(which(loaded_metadata[metadata_category,] == group_1))
      group_2_samples <- names(which(loaded_metadata[metadata_category,] == group_2))
      group_interaction_score_overview_pvalues[interaction,comparison] <- t.test(interaction_score_overview[interaction,group_1_samples],interaction_score_overview[interaction,group_2_samples])$p.value
    }
  }
  # Store interaction overview, correlations, group overview, and group p values
  k_interaction_score_overview[[as.character(interaction_k)]] <- interaction_score_overview
  k_group_interaction_score_overview[[as.character(interaction_k)]] <- group_interaction_score_overview
  k_group_interaction_score_overview_pvalues[[as.character(interaction_k)]] <- group_interaction_score_overview_pvalues
}
# Populate differential networks with -log10(p values) for visualization 
k_differential_pvalue_color <- k_differential_interaction_color
k_differential_pvalue_scores <- list()
k_differential_pvalue_network <- list()
for (interaction_k in c(6,18,36,60,90))
{
  # Load variables for current k parameter
  group_interaction_score_overview_pvalues <- k_group_interaction_score_overview_pvalues[[as.character(interaction_k)]]
  # For discrete groups, populate networks with appropriate -log10(p values)
  differential_pvalue_scores <- list()
  differential_pvalue_network <- list()
  for (comparison in discrete_pairs)
  {
    # Define variable for storing values
    differential_pvalue_scores[[comparison]] <- matrix(nrow = nrow(cell_pairs), ncol=3)
    colnames(differential_pvalue_scores[[comparison]]) <- c("f1","f2","score") 
    differential_pvalue_scores[[comparison]][,"f1"] <- cell_pairs[,1]
    differential_pvalue_scores[[comparison]][,"f2"] <- cell_pairs[,2]
    rownames(differential_pvalue_scores[[comparison]]) <- paste(differential_pvalue_scores[[comparison]][,"f1"], differential_pvalue_scores[[comparison]][,"f2"], sep = " - ")
    # Populate values
    differential_pvalue_scores[[comparison]][,"score"] <- -log10(group_interaction_score_overview_pvalues[rownames(differential_pvalue_scores[[comparison]]),comparison])
    # Generate differential interaction score matrix
    differential_pvalue_score_matrix <- AdjScoreHeatmapMatrix(differential_pvalue_scores[[comparison]])
    class(differential_pvalue_score_matrix) <- "numeric"
    # Generate differential interaction network
    differential_pvalue_network[[comparison]] <- graph_from_adjacency_matrix(differential_pvalue_score_matrix, weighted = TRUE)
  }
  # Store differential interaction scores, networks, and colors
  k_differential_pvalue_scores[[as.character(interaction_k)]] <- differential_pvalue_scores
  k_differential_pvalue_network[[as.character(interaction_k)]] <- differential_pvalue_network
  rm(differential_pvalue_score_matrix)
}
rm(differential_pvalue_scores,differential_pvalue_network)
# Plot differential interaction network
interaction_k <- 6 # Specify k of interest
comparison_of_interest <- "D14.Group vs. D35.Group" # Specify comparison of interest
# Populate variables for plotting
differential_pvalue_network <- k_differential_pvalue_network[[as.character(interaction_k)]]
differential_pvalue_color <- k_differential_pvalue_color[[as.character(interaction_k)]]
# Plot interaction network for specific comparison
ggraph(differential_pvalue_network[[comparison_of_interest]], layout = "linear", circular = TRUE) + geom_edge_link(aes(edge_width=E(differential_pvalue_network[[comparison_of_interest]])$weight, edge_colour = as.factor(as.vector(differential_pvalue_color[[comparison_of_interest]][differential_pvalue_color[[comparison_of_interest]] != 0])))) + geom_node_point( color="black", size=8) + geom_node_text( aes(label=name), repel = TRUE, size=8, color="#681984") + theme_void() + theme(legend.position="none", plot.margin=unit(rep(1,4), "cm")) + scale_edge_colour_manual(values = c("#317cca48","#C43d3d48"))

