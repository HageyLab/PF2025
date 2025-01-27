## Adapted from Mascharak & Guo et al. Cell Reports Medicine 2023.

## Initialization
library(STvEA)
library(AdjacencyScore)
library(RANN)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(patchwork)
library(monocle3)
library(ggplot2)
library(gridExtra)
library(stringr)
library(flowCore)
library(batchelor)
library(igraph)
library(ggraph)
library(mixOmics)
library(Cairo)
source("Additional CODEX Functions.R")


# Prior to this point, export the .fcs files for each imaged CODEX region and concatenate into a single file (all_regions.fcs) using flowJo, Excel, or equivalent software.
# Additionally, define stain gating in the native CODEX software (multiplex analysis viewer; MAV) and export as a single file (gate.csv)


## Pre-Processing Raw Data in STvEA (Govek et al. Science Advances 2019.)
# Read in channel names from .fcs file
codex_names <- ReadNamesFCS("all_regions.fcs")[["channels"]]
# Store names of blank, non-protein, and DAPI "channels"
is_blank <- grepl("Blank",codex_names)
is_protein <- !(codex_names %in% c("cell_id","region","tile_num","x","y","z","x_tile","y_tile","size","tsne_x","tsne_y","homogeneity"))
is_DAPI <- grepl("DAPI",codex_names)
# Read in dataset with information defined above
stvea_object <- ReadDataFCS("all_regions.fcs", is_protein, is_blank, protein_names = codex_names[is_protein])
# Read in MAV-defined gating
gate <- read.delim(file = "gate.csv", header = TRUE, sep = ",")
gate <- as.logical(unlist(gate[,2]))
# Filter out stained debris using method from STvEA paper and MAV-defined gating
stvea_object <- FilterCODEX(stvea_object, inclusion = gate)
# Normalize staining intensities to Gaussian distribution
stvea_object <- CleanCODEX(stvea_object, model="gaussian", num_cores=4)


## Batch Correction in mnnCorrect (Haghverdi et al. Nature Biotechnology 2018.)
# Set batch correction parameters
b <- 200
# Perform batch correction between image regions using mutual nearest neighbors
temp_object <- mnnCorrect(t(stvea_object@codex_clean), batch = stvea_object@codex_region, k = b, cos.norm.in = FALSE, cos.norm.out = FALSE, subset.row = !grepl("DAPI",codex_names[is_protein]) & !grepl("Blank",codex_names[is_protein]), correct.all = TRUE)
stvea_object@codex_clean <- t(temp_object@assays@data@listData$corrected)
rm(temp_object)
gc()


## Generating UMAP Manifold
# Generate counts matrix for testing dataset
seurat_matrix <- stvea_object@codex_clean[,!grepl("Blank",colnames(stvea_object@codex_clean))]
rownames(seurat_matrix) <- paste(stvea_object@codex_region, stvea_object@codex_cell_ID, sep = "_")
seurat_matrix <- t(seurat_matrix)
# Create Seurat object for testing dataset
# Set UMAP parameters of interest
m <- "cosine"; i <- 70; j <- 0.1; k <- 1; l <- 100
# Create Seurat object and run minimum required preprocessing
seurat_object <- CreateSeuratObject(counts = seurat_matrix)
rm(seurat_matrix)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- RunPCA(seurat_object, features = rownames(seurat_object), npcs = 5)
# UMAP processing (timed)
umap_start_time <- Sys.time()
seurat_object <- RunUMAP(seurat_object, umap.method = "umap-learn", return.model = TRUE, dims = 1:5, metric = m, n.neighbors = i, min.dist = j, repulsion.strength = k, negative.sample.rate = l)
umap_end_time <- Sys.time()
# Generate UMAP clusters
seurat_object <- FindNeighbors(seurat_object, dims = 1:5)
seurat_object <- FindClusters(seurat_object, resolution = 0.5) 
# Plot UMAP manifold
DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + NoLegend()
gc()


## Transferring Information Back to STvEA
# Populate STvEA object with Seurat UMAP, clusters, etc.
stvea_object@codex_emb <- data.frame(row.names = rownames(stvea_object@codex_protein))
stvea_object@codex_emb$V1 <- seurat_object[["umap"]]@cell.embeddings[,1]
stvea_object@codex_emb$V2 <- seurat_object[["umap"]]@cell.embeddings[,2]
stvea_object@codex_clusters <- seurat_object$seurat_clusters
# Plot manifold in STvEA
PlotClusterCODEXemb(stvea_object)


## Cell Type Annotations
# Import cell type annotations into STvEA
loaded_annotations <- read.delim(file = "cluster_annotations.csv", header = TRUE, row.names=1, sep = ",")
loaded_annotations <- t(loaded_annotations)
stvea_object@codex_clusters <- as.factor(loaded_annotations[,as.character(stvea_object@codex_clusters)])
# Add cell type annotations to Seurat
seurat_object$annotated_cell_types <- as.character(stvea_object@codex_clusters)


## Interaction Analysis
# Generate cell pairings
cell_types <- sort(as.character(unique(stvea_object@codex_clusters)))
individual_cell_pairs <- t(combn(cell_types,2))
individual_cell_pairs <- rbind(individual_cell_pairs,cbind(individual_cell_pairs[,2],individual_cell_pairs[,1]))
for (cell in cell_types) 
{
  individual_cell_pairs <- rbind(individual_cell_pairs, c(cell,cell))
}
# Set k parameter
interaction_k <- 20
# Generate interaction scores
individual_interaction_scores <- list()
for (region in as.character(unique(stvea_object@codex_region)))
{
  # Generate adjacency matrix based on cell locations
  adjacency_matrix <- knn_graph(stvea_object@codex_spatial[stvea_object@codex_region == region,], k=interaction_k)
  # Generate "feature" matrix of assigned cell types
  cell_matrix <- t(sapply(cell_types, function(x) (stvea_object@codex_clusters[stvea_object@codex_region == region]==x)*1))
  row.names(cell_matrix) <- cell_types
  # Calculate interaction scores as dot-product of paired cell feature vectors with adjacency matrix (f1*j*f2)
  individual_interaction_scores[[region]] <- CalculateInteractionScores(adjacency_matrix = adjacency_matrix, cell_matrix = cell_matrix, cell_pairs = individual_cell_pairs, calculate_individual_scores = TRUE, normalize = FALSE, cell_labels = colnames(seurat_object)[stvea_object@codex_region == region])
}
rm(adjacency_matrix)
# Store interaction scores 
cellular_interaction_scores <- NULL
for (sample in names(individual_interaction_scores))
{
  if (is.null(cellular_interaction_scores))
  {
    cellular_interaction_scores <- individual_interaction_scores[[sample]]
  } else
  {
    cellular_interaction_scores <- cbind(cellular_interaction_scores, individual_interaction_scores[[sample]])
  }
}
cellular_interaction_scores <- cellular_interaction_scores[,!(colnames(cellular_interaction_scores) %in% c("f1","f2"))]


