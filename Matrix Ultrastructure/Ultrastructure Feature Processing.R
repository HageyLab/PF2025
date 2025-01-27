## Adapted from Mascharak & Guo et al. Cell Reports Medicine 2023.

## Initialization
library(dplyr)
library(patchwork)
library(DDRTree)
library(ggplot2)
library(monocle3)
library(gridExtra)


## Loading in MATLAB-Processed Ultrastructure Data
# Load data (.csv should have one row per image and one column per matrix parameter and UMAP dimension, plus a header column and a group column)
loaded_counts <- read.delim(file = "Quantified.csv", header = TRUE, row.names=1, sep = ",")
loaded_counts <- t(loaded_counts)
# Store group key
loaded_key <- loaded_counts["Group",]
# Store UMAP coordinates
loaded_UMAP <- matrix(nrow = 2, ncol = ncol(loaded_counts))
loaded_UMAP[1,] <- as.numeric(loaded_counts[rownames(loaded_counts) == "UMAP.X",]) 
loaded_UMAP[2,] <- as.numeric(loaded_counts[rownames(loaded_counts) == "UMAP.Y",]) 
# Store parameter data
loaded_counts <- loaded_counts[rownames(loaded_counts)!="Group" & !grepl("UMAP",rownames(loaded_counts)),]


## Trajectory Analysis (DDRTree via Monocle Wrappers)
# Initialize CDS object
us_object <- new_cell_data_set(loaded_counts)
# Store group assignments
colData(us_object)$group <- as.character(loaded_key)
# Store parameters
rowData(us_object)$gene_name <- rownames(us_object)
rowData(us_object)$gene_short_name <- rowData(us_object)$gene_name
# Run UMAP function to generate dummy fields for import of MATLAB-processed coordinates
us_object <- preprocess_cds(cds = us_object, method = "PCA", num_dim=5)
us_object <- reduce_dimension(cds = us_object, reduction_method = "UMAP", preprocess_method="PCA")
reducedDim(us_object, "UMAP")[,1] <- loaded_UMAP[1,]
reducedDim(us_object, "UMAP")[,2] <- loaded_UMAP[2,]
# Build trajectories (DDRTree via Monocle)
# Optional: resolution and minimal_branch_len may be tuned to learn a graph with a maximum degree of 2
# This helps eliminate excess branches, but does not significantly impact assigned pseudotime values
us_object <- cluster_cells(cds = us_object, reduction_method = "UMAP", resolution = 3e-4) 
us_object <- learn_graph(us_object, use_partition = TRUE, learn_graph_control = list(minimal_branch_len = 15))
us_object <- order_cells(us_object, reduction_method = "UMAP")
# Plot pseudotime trajectories
plot_cells(cds = us_object, color_cells_by = "pseudotime", cell_size = 2, show_trajectory_graph = TRUE, trajectory_graph_segment_size = 1.5, label_branch_points = FALSE, graph_label_size = 4, group_label_size = 0) + guides(colour=guide_colourbar(title='Pseudotime', ticks = FALSE))

