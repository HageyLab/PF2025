## Adapted from  Mascharak & Guo et al. Cell Reports Medicine 2023.

## Initialization
library(magick)
library(RANN)
library(monocle3)
library(AdjacencyScore)
library(igraph)
library(ggraph)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(gridExtra)
source("Additional Functions.R")


# Prior to this point, place .tif files for each Visium-associated WSI in the
# current working directory, create a subfolder named "Tiles" for storage of 
# programamatically generated tiles, and transfer each SampleID_tissue_positions.csv file 
# from the aligned Visium data to the current working directory. Additionally, load 
# a reference CDS (us_object) containing  the originally trained pseudotime model, 
# as well as all objects generated by "Cell Spatial Interaction Analysis.R".


## WSI Tiling
# Specify desired tile width and height (choose pixel dimensions that represent
# desired physical length scale, e.g. width/height in microns used in original analysis)
tile_width <- 389
tile_height <- 293
# Import test image and store metadata
for (sample_name in c("A1","A2","B1","B2","C1","C2","D1","D2"))
{
  whole_image <- image_read(paste0(getwd(),"/",sample_name,".tif"))
  whole_image_info <- image_info(whole_image)
  # Store and export tiles
  for (x in seq(0,whole_image_info$width%/%tile_width*tile_width,tile_width))
  {
    for (y in seq(0,whole_image_info$height%/%tile_height*tile_height,tile_height))
    {
      # Crop and export tile
      tiled_image <- image_crop(whole_image, sprintf(paste0(tile_width,"x",tile_height,"+%d+%d"), x, y))
      image_write(tiled_image, path = paste0(getwd(),"/Tiles/",sample_name,"_",x,"_",y,".tif"), format = "tif")
    }
  }
}


# At this point, remove any tiles smaller than the desired dimensions (i.e. edge tiles) and run
# ultrastructure analysis in MATLAB on all remaining tiles. Export quantified.csv.


## Importing Original Pixel-Level Coordinates from Visium
# Load and format original coordinates from full-resolution images
original_coordinates <- list()
for (sample in c("A1","A2","B1","B2","C1","C2","D1","D2"))
{
  # Read in coordinates from Cellranger .csv file (limit to spots included in capture)
  sample_coordinates <- read.delim(file = paste0(getwd(),"/",sample,"_tissue_positions.csv"), header = TRUE, row.names=1, sep = ",")
  sample_coordinates <- sample_coordinates[sample_coordinates[,"in_tissue"] == 1,]
  sample_coordinates <- sample_coordinates[,colnames(sample_coordinates) %in% c("pxl_row_in_fullres","pxl_col_in_fullres")]
  # Change to standardized coordinate system
  colnames(sample_coordinates)[colnames(sample_coordinates) == "pxl_row_in_fullres"] <- "y_centroid"
  colnames(sample_coordinates)[colnames(sample_coordinates) == "pxl_col_in_fullres"] <- "x_centroid"
  # Format spot ID's 
  rownames(sample_coordinates) <- paste0(sample,"_",rownames(sample_coordinates)) 
  # Store coordinates
  original_coordinates[[sample]] <- sample_coordinates
  rm(sample_coordinates)
}


## Mapping Ultrastructure Data in Monocle3 and RANN
# Specify tile width and height
tile_width <- 389
tile_height <- 293
# Load ultrastructure data (.csv should have one row per image and one column per matrix parameter, plus headers)
test_counts <- read.delim(file = "quantified.csv", header = TRUE, row.names=1, sep = ",")
test_counts[is.na(test_counts)] <- 0
tile_data <- test_counts
# Add tile metadata
tile_data$sample <- NA
tile_data$lower_x <- NA
tile_data$lower_y <- NA
tile_data$upper_x <- NA
tile_data$upper_y <- NA
for (tile in rownames(tile_data))
{
  split_name <- strsplit(tile, split = "_")[[1]]
  tile_data[tile,"sample"] <- split_name[1]
  tile_data[tile,"lower_x"] <- as.numeric(split_name[2])
  tile_data[tile,"lower_y"] <- as.numeric(split_name[3])
  tile_data[tile,"upper_x"] <- tile_data[tile,"lower_x"] + tile_width
  tile_data[tile,"upper_y"] <- tile_data[tile,"lower_y"] + tile_height
}
# Find number of spots that fall within each tile
tile_data$n_spots <- NA
for (tile in rownames(tile_data))
{
  # Load original pixel coordinates for sample
  sample_coordinates <- original_coordinates[[tile_data[tile,"sample"]]]
  # Find spots that fall within pixel coordinates of tile
  spots <- rownames(sample_coordinates[sample_coordinates$x_centroid > tile_data[tile,"lower_x"] & sample_coordinates$x_centroid <= tile_data[tile,"upper_x"] & sample_coordinates$y_centroid > tile_data[tile,"lower_y"] & sample_coordinates$y_centroid <= tile_data[tile,"upper_y"],])
  # Store number of spots in tile
  tile_data[tile,"n_spots"] <- length(spots)
  rm(sample_coordinates,spots)
}
# Subset tiles to only include those with Visium spots
test_counts <- test_counts[tile_data$n_spots >= 1,]
# Format counts for downstream import
test_counts <- t(test_counts)
# Convert single color (e.g. HE) parameters to equivalent double color (e.g. Picrosirius Red) parameters
parameter_key <- read.delim(file = "single_to_double_color_parameter_key.csv", header = TRUE, row.names=1, sep = ",")
rownames(test_counts) <- parameter_key[rownames(test_counts),"double_color_parameter"]
# Create objects for nearest neighbor indexing and new query data
nn_object <- us_object
test_object <- new_cell_data_set(as(test_counts,"dgCMatrix"))
rowData(test_object)$gene_name <- rownames(test_object)
rowData(test_object)$gene_short_name <- rowData(test_object)$gene_name
# Limit original and new object to shared matrix parameters
shared_parameters <- intersect(rownames(test_counts),rownames(nn_object))
test_object <- test_object[shared_parameters,]
nn_object <- nn_object[shared_parameters,]
# Generate nearest neighbor index for original dataset 
nn_object <- preprocess_cds(cds = nn_object, method = "PCA", num_dim=100)
nn_object <- reduce_dimension(cds = nn_object, reduction_method = "UMAP", preprocess_method="PCA", build_nn_index = TRUE)
nn_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- us_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
# Extract model from original dataset and load into query dataset
save_transform_models(nn_object, 'reference_model')
test_object <- load_transform_models(test_object, "reference_model")
# Apply original transform model to query dataset
test_object <- preprocess_transform(test_object)
test_object <- reduce_dimension_transform(test_object)
# Map pseudotime values to new images by nearest neighbors
nn_indices <- as.vector(nn2(data = reducedDim(nn_object, "PCA"), query = reducedDim(test_object, "PCA"), k = 1)$nn.idx)
test_object$mapped_pseudotime <- nn_object@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][as.numeric(nn_indices)]
mapped_pseudotime <- test_object$mapped_pseudotime
# Export matrix parameters and mapped pseudotime values
mapped_values <- rbind(test_counts,mapped_pseudotime)
rownames(mapped_values)[nrow(mapped_values)] <- "Mapped.Pseudotime"
write.csv(t(mapped_values),"quantified_mapped.csv")
# Recreate tile dataframe and associated variables based on mapped values
tile_data <- read.delim(file = "quantified_mapped.csv", header = TRUE, row.names=1, sep = ",")
tile_data$sample <- NA
tile_data$lower_x <- NA
tile_data$lower_y <- NA
tile_data$upper_x <- NA
tile_data$upper_y <- NA
for (tile in rownames(tile_data))
{
  split_name <- strsplit(tile, split = "_")[[1]]
  tile_data[tile,"sample"] <- split_name[1]
  tile_data[tile,"lower_x"] <- as.numeric(split_name[2])
  tile_data[tile,"lower_y"] <- as.numeric(split_name[3])
  tile_data[tile,"upper_x"] <- tile_data[tile,"lower_x"] + tile_width
  tile_data[tile,"upper_y"] <- tile_data[tile,"lower_y"] + tile_height
}
tile_data$n_spots <- NA
tile_spots <- list()
for (tile in rownames(tile_data))
{
  # Load original pixel coordinates for sample
  sample_coordinates <- original_coordinates[[tile_data[tile,"sample"]]]
  # Find spots that fall within pixel coordinates of tile 
  spots <- rownames(sample_coordinates[sample_coordinates$x_centroid > tile_data[tile,"lower_x"] & sample_coordinates$x_centroid <= tile_data[tile,"upper_x"] & sample_coordinates$y_centroid > tile_data[tile,"lower_y"] & sample_coordinates$y_centroid <= tile_data[tile,"upper_y"],])
  spots <- spots[spots %in% colnames(vsm)]
  # Store spot ID's 
  tile_spots[[tile]] <- spots 
  # Store number of spots in tile
  tile_data[tile,"n_spots"] <- length(tile_spots[[tile]])
  rm(sample_coordinates,spots)
}
# Store matrix parameters (with pseudotime listed first) for downstream utilization
matrix_parameters <- colnames(tile_data)
matrix_parameters <- c("Mapped.Pseudotime",matrix_parameters[matrix_parameters != "Mapped.Pseudotime"])


## Tile-Level Visium Analysis
# Initialize variables for storage of cell probabilities
relevant_cell_types <- rownames(vsm@assays$predictions@data)[rownames(vsm@assays$predictions@data) != "max"]
for (cell_type in relevant_cell_types)
{
  tile_data[[cell_type]] <- NA
}
# Calculate normalized cell probabilities in each tile
for (tile in rownames(tile_data))
{
  # Load spot ID's within tile
  spots <- tile_spots[[tile]]
  # Average cell probabilities for spots within tile
  if (length(spots) == 1)
  {
    tile_data[tile,relevant_cell_types] <-vsm@assays$predictions@data[relevant_cell_types,spots]
  } else
  {
    tile_data[tile,relevant_cell_types] <- rowSums(vsm@assays$predictions@data[relevant_cell_types,spots])/length(spots)
  }
}
# Initialize variables for storage of cell interactions
relevant_interactions <- rownames(empty_matrix)
for (interaction in relevant_interactions)
{
  tile_data[[interaction]] <- NA
}
interaction_k <- 6
# Calculate normalized interaction scores in each tile
for (tile in rownames(tile_data))
{
  # If too few spots, populate with empty matrix
  if (tile_data[tile,"n_spots"] < interaction_k)
  {
    tile_int_scores <- empty_matrix
  } else
  {
    # Generate adjacency matrix based on cell locations
    spots <- tile_spots[[tile]]
    sample <- names(which(loaded_metadata["Sample.Name",] == tile_data[tile,"sample"]))
    adjacency_matrix <- knn_graph(cell_locations[[sample]][spots,], k = interaction_k)
    # Calculate interaction scores as dot-product of paired cell feature vectors with adjacency matrix (f1*j*f2)
    tile_int_scores <- CalculateInteractionScores(adjacency_matrix = adjacency_matrix, cell_matrix = t(cell_probability_matrix[[sample]][spots,]), cell_pairs = cell_pairs, calculate_individual_scores = FALSE, normalize = TRUE)
    rm(adjacency_matrix)
  }
  # Store tile-level interaction scores
  tile_data[tile,relevant_interactions] <- tile_int_scores[relevant_interactions,"score"]
}
# Initialize variables for storage of gene expression
norm_matrix <- vsm@assays$SCT@data
for (gene in rownames(norm_matrix))
{
  tile_data[[gene]] <- NA
}
# Calculate normalized gene expression in each tile
for (tile in rownames(tile_data))
{
  # Load spot ID's within tile
  spots <- tile_spots[[tile]]
  # Average gene expression for spots within tile
  if (length(spots) == 1)
  {
    tile_data[tile,rownames(norm_matrix)] <- norm_matrix[,spots]
  } else
  {
    tile_data[tile,rownames(norm_matrix)] <- rowSums(norm_matrix[,spots])/length(spots)
  }
}
# Export all tile-level data for external analysis
write.csv(tile_data,"matrix_visium_tile_data.csv")
# Populate matrix parameters in Visium object for visualization
spots_matrix <- matrix(nrow = length(matrix_parameters), ncol = ncol(vsm))
rownames(spots_matrix) <- matrix_parameters
colnames(spots_matrix) <- colnames(vsm)
for (parameter in matrix_parameters)
{
  for (tile in names(tile_spots))
  {
    spots_matrix[parameter,tile_spots[[tile]]] <- tile_data[tile,parameter]
  }
  vsm@meta.data[[parameter]] <- spots_matrix[parameter,]
}
rm(spots_matrix)
# Plot matrix pseudotime on Visium
SpatialFeaturePlot(vsm, features = "Mapped.Pseudotime", images = "slice1_B1")


