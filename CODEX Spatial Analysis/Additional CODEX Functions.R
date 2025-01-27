#' Create version of STvEA class (Govek et al. Science Advances 2019.)
#' that includes slots for cell ID and region
#' Class to hold all data matrices needed for mapping CODEX
#' and following analysis
#'
STvEA.data <- setClass(
  Class = 'STvEA.data',
  slots = c(
    cite_mRNA = 'ANY', # cite_cells x genes
    cite_mRNA_norm = 'ANY', # cite cells x genes - log(1+TPM) or seurat
    cite_latent = 'ANY', # cite cells x low
    cite_emb = 'ANY', # cite cells x 2
    hdbscan_param_scan = 'list', # output of ParameterScan
    cite_clusters = 'vector',
    cite_protein = 'ANY', # cite cells x proteins
    cite_clean = 'ANY', # cite cells x proteins
    cite_norm = 'ANY', # cite cells x proteins

    codex_protein = 'ANY', # codex cells x proteins
    codex_size = 'numeric',
    codex_cell_ID = 'vector',
    codex_region = 'vector',
    codex_blanks = 'ANY', # codex cells x blank channels
    codex_spatial = 'ANY', # codex cells x 3
    codex_emb = 'ANY', # codex cells x 2
    codex_knn = 'ANY', # codex cells x k
    codex_clusters = 'vector',
    codex_clean = 'ANY', # codex cells x proteins

    corrected_codex = 'ANY', # codex cells x proteins
    transfer_matrix = 'dgCMatrix', # codex cells x cite cells
    # k CODEX neighbors for each CITE-seq cell
    codex_mRNA = 'ANY' # codex cells x genes
  )
)

#' Create version of SetDataCODEX function (Govek et al. Science Advances 2019.)
#' that includes cell ID and region information
#' Set CODEX data in STvEA.data object
#'
#' @param codex_protein Segmented and spillover corrected CODEX protein expression (cell x protein)
#' @param codex_blanks Segmented and spillover corrected CODEX blank channel expression (cell x blank channel)
#' @param codex_size Size of each CODEX cell (vector)
#' @param codex_spatial (optional) xyz coordinates of each CODEX cell (cell x 3 coordinates)
#' @param codex_cell_ID (optional) ID number of each CODEX cell (vector)
#' @param codex_region (optional) region number of each CODEX cell (vector)
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @export
#'
SetDataCODEX <- function(codex_protein,
                               codex_blanks,
                               codex_size,
                               codex_spatial = NULL,
                               codex_cell_ID = NULL,
                               codex_region = NULL,
                               stvea_object = NULL) {
  if (is.null(stvea_object)) {
    stvea_object <- new(
      Class = "STvEA.data",
      codex_protein = codex_protein,
      codex_blanks = codex_blanks,
      codex_size = codex_size,
      codex_spatial = codex_spatial,
      codex_cell_ID = codex_cell_ID,
      codex_region = codex_region
    )
  } else {
    stvea_object@codex_protein <- codex_protein
    stvea_object@codex_blanks <- codex_blanks
    stvea_object@codex_size <- codex_size
    stvea_object@codex_spatial <- codex_spatial
    stvea_object@codex_cell_ID <- codex_cell_ID
    stvea_object@codex_region <- codex_region
  }
  return(stvea_object)
}

#' Create version of ReadDataFCS function that reads in cell ID and
#' region information from the .fcs file
#' Read in FCS file
#'
#' @param path_to_fcs readable path to FCS file
#' @param is_protein boolean vector indicating which columns in FCS
#' should be kept as protein expression for downstream analyses
#' @param is_blank boolean vector indicating which columns in FCS
#' are blank for filtering purposes
#' @param protein_names vector of names for each protein channel column -
#' if null, will keep column names from FCS
#' @param tile_relative whether x,y,z coordinates are relative to x_tile,y_tile (TRUE)
#' or to corner of image (FALSE).
#'       If TRUE, will compute absolute coordinates
#'       If FALSE, will just take given x,y,z coordinates
#' @param num_tile_x number of tiles in the x direction, only used if tile_relative = TRUE
#' and the tile information in the FCS file is saved as tile_nr instead of x_tile,y_tile
#' @param x_pix_size x dimension of pixel coordinates in nm. Set to 1 to keep pixel dimensions.
#' @param y_pix_size y dimension of pixel coordinates in nm. Set to 1 to keep pixel dimensions.
#' @param z_pix_size z dimension of pixel coordinates in nm. Set to 1 to keep pixel dimensions.
#' @param stvea_object (optional) Pre-existing STvEA.data object to load data into.
#' If not provided, a new object is created.
#'
#' @return STvEA.data class object
#'
#' @importFrom flowCore read.FCS
#' @importFrom stringr str_match
#'
#' @export
#'
ReadDataFCS <- function(path_to_fcs,
                              is_protein,
                              is_blank,
                              protein_names = NULL,
                              tile_relative = FALSE,
                              num_tile_x = 9,
                              x_pix_size = 188,
                              y_pix_size = 188,
                              z_pix_size = 900,
                              stvea_object = NULL) {
  expr_mat <- read.FCS(path_to_fcs,transformation=FALSE, truncate_max_range=FALSE)@exprs

  # Get blank columns
  if (sum(is_blank) == 0) {
    warning("No columns are labelled as blank, setting blank channel data to all 0s")
    codex_blanks <- as.matrix(rep(0,nrow(expr_mat)), ncol=1)
  } else {
    codex_blanks <- expr_mat[,is_blank,drop=FALSE]
  }

  # Get protein matrix
  if (!is.null(protein_names) && sum(is_protein) != length(protein_names)) {
    stop("Length of protein_names must be the same as TRUE values in is_protein")
  }
  codex_protein <- expr_mat[,is_protein,drop=FALSE]
  if (!is.null(protein_names)) {
    colnames(codex_protein) <- protein_names
  }

  # Match columns of the format "X.X" and convert them to just "x"
  reg_match <- str_match(colnames(expr_mat), "(.+)\\.\\1")
  match <- !is.na(reg_match[,1]) # these columns have the form "X.X"
  colnames(expr_mat)[match] <- reg_match[match,2] # convert to just "X"
  colnames(expr_mat) <- sapply(colnames(expr_mat), tolower) # lowercase "x"

  # Get spatial information
  if (!all(c("x","y","z") %in% colnames(expr_mat))) {
    warning("Cannot find x,y,z coordinates in FCS. Continuing without spatial information")
    codex_spatial_nm <- NULL
  } else {
    if (tile_relative) {
      if ("tile_nr" %in% colnames(expr_mat)) {
        x <- floor((expr_mat[,"tile_nr"]-1)/num_tile_x) * max(expr_mat[,"x"]) + expr_mat[,"x"]
        y <- ((expr_mat[,"tile_nr"] - 1) %% num_tile_x) * max(expr_mat[,"y"]) + expr_mat[,"y"]
      } else if ("tile_num" %in% colnames(expr_mat)) {
        x <- floor((expr_mat[,"tile_num"]-1)/num_tile_x) * max(expr_mat[,"x"]) + expr_mat[,"x"]
        y <- ((expr_mat[,"tile_num"] - 1) %% num_tile_x) * max(expr_mat[,"y"]) + expr_mat[,"y"]
      } else if ("x_tile" %in% colnames(expr_mat) && "y_tile" %in% colnames(expr_mat)) {
        x <- (expr_mat[,"x_tile"] - 1) * max(expr_mat[,"x"]) + expr_mat[,"x"]
        y <- (expr_mat[,"y_tile"] - 1) * max(expr_mat[,"y"]) + expr_mat[,"y"]
      } else {
        stop("Cannot find tile information in FCS to compute absolute x,y,z from relative")
      }
      z <- expr_mat[,"z"]
      codex_spatial <- cbind(x,y,z)
    } else {
      codex_spatial <- expr_mat[,c("x","y","z")]
    }
    codex_spatial_nm <- cbind(x=codex_spatial[,"x"]*x_pix_size,
                              y=codex_spatial[,"y"]*y_pix_size,
                              z=codex_spatial[,"z"]*z_pix_size)
  }

  # Get size information
  if ("size" %in% colnames(expr_mat)) {
    codex_size <- expr_mat[,"size"]
  } else {
    warning("Cannot find size information in FCS, setting all cell size to 0")
    codex_size <- rep(0, nrow(expr_mat))
  }

  if ("cell_id" %in% colnames(expr_mat))
  {
    codex_cell_ID <- expr_mat[,"cell_id"]
  } else {
    warning("Cannot find cell IDs in FCS, setting cell IDs to sequential numbers")
    codex_cell_ID <- seq(length.out = nrow(expr_mat))
  }
  if ("region" %in% colnames(expr_mat))
  {
    codex_region <- expr_mat[,"region"]
  } else {
    warning("Cannot find regions in FCS, setting all regions to 0")
    codex_region <- rep(0, nrow(expr_mat))
  }

  stvea_object <- SetDataCODEX(codex_protein = codex_protein,
                                     codex_blanks = codex_blanks,
                                     codex_size = codex_size,
                                     codex_spatial = codex_spatial_nm,
                                     codex_cell_ID = codex_cell_ID,
                                     codex_region = codex_region,
                                     stvea_object = stvea_object)

  return(stvea_object)
}

#' Subsetting of CODEX object to specified region IDs
#'
#' @param stvea_object Pre-existing STvEA.data object
#' @param region_IDs Charcater vector of region IDs (character) to keep
#'
#' @export
#'
SubsetCODEXRegions <- function(stvea_object,
                         region_IDs) {
  if (is.null(stvea_object))
  {
    stop("Please supply an existing STvEA object.", call. =FALSE)
  }
  if (is.null(region_IDs) | !is.character(region_IDs) | !is.vector(region_IDs))
  {
    stop("Please supply a character vector of region ID's to keep.", call. =FALSE)
  }
  cells_to_keep <- stvea_object@codex_region %in% region_IDs
  if (!is.null(stvea_object@codex_protein))
  {
    stvea_object@codex_protein <- stvea_object@codex_protein[cells_to_keep,]
  }
  if (!is.null(stvea_object@codex_size))
  {
    stvea_object@codex_size <- stvea_object@codex_size[cells_to_keep]
  }
  if (!is.null(stvea_object@codex_cell_ID))
  {
    stvea_object@codex_cell_ID <- stvea_object@codex_cell_ID[cells_to_keep]
  }
  if (!is.null(stvea_object@codex_region))
  {
    stvea_object@codex_region <- stvea_object@codex_region[cells_to_keep]
  }
  if (!is.null(stvea_object@codex_blanks))
  {
    stvea_object@codex_blanks <- stvea_object@codex_blanks[cells_to_keep,]
  }
  if (!is.null(stvea_object@codex_spatial))
  {
    stvea_object@codex_spatial <- stvea_object@codex_spatial[cells_to_keep,]
  }
  if (!is.null(stvea_object@codex_emb))
  {
    stvea_object@codex_emb <- stvea_object@codex_emb[cells_to_keep,]
  }
  if (!is.null(stvea_object@codex_knn))
  {
    stvea_object@codex_knn <- stvea_object@codex_knn[cells_to_keep,]
  }
  if (!is.null(stvea_object@codex_clusters))
  {
    stvea_object@codex_clusters <- stvea_object@codex_clusters[cells_to_keep]
  }
  if (!is.null(stvea_object@codex_clean))
  {
    stvea_object@codex_clean <- stvea_object@codex_clean[cells_to_keep,]
  }
  return(stvea_object)
}

#' Subsamples each CODEX region to a either specified fraction of the original number
#' of cells in each region, or to a fixed number of cells per region
#'
#' @param stvea_object STvEA.data class with CODEX data
#' @param use_fraction a boolean, with TRUE indicating usage of fraction and FALSE indicating usage of fixed number
#' @param random_sampling either the fraction of cells to select from each region or a fixed number of cells to select from each region, depending on use_fraction 
#'
#' @export
#'
SubsampleCODEX <- function(stvea_object,
                              use_fraction = TRUE,
                              random_sampling = 0.1) {
  # Check that stvea_object parameter is the correct class
  if (class(stvea_object) != "STvEA.data")
  {
    stop("stvea_object parameter must be of class STvEA.data", call. =FALSE)
  }
  # Check that use_fraction parameter is a boolean
  if (!is.logical(use_fraction) | !is.vector(use_fraction) | length(use_fraction) != 1)
  {
    stop("use_fraction parameter must be a single boolean.", call. =FALSE)
  }
  # Check that random_sampling parameter is set correctly in call to function
  if (!is.numeric(random_sampling) | !is.vector(random_sampling) | length(random_sampling) != 1)
  {
    stop("random_sampling parameter must be a single number.", call. =FALSE)
  }
  
  # If specified by use_fraction = TRUE, calculate per-region number of cells as a fraction of original number
  if (use_fraction)
  {
    if (random_sampling >= 1)
    {
      stop("If use_fraction = TRUE, random_sampling must be < 1.", call. =FALSE)
    }
    # Determine number of points to keep in each region
    region_sizes <- vector(mode = "integer", length = length(unique(stvea_object@codex_region)))
    names(region_sizes) <- unique(stvea_object@codex_region)
    for (region in names(region_sizes))
    {
      individual_region_size <- sum(stvea_object@codex_region == region)
      region_sizes[region] <- round(individual_region_size*random_sampling)
      if (region_sizes[region]==0)
      {
        region_sizes[region] <- 1
      }
    }
  }
  
  # If specified by use_fraction = FALSE, specify subsampling using a fixed number of cells
  if (!use_fraction)
  {
    random_sampling <- round(random_sampling)
    # Determine number of points to keep in each region
    region_sizes <- vector(mode = "integer", length = length(unique(stvea_object@codex_region)))
    names(region_sizes) <- unique(stvea_object@codex_region)
    for (region in names(region_sizes))
    {
      if (random_sampling > sum(stvea_object@codex_region == region))
      {
        stop("If use_fraction = FALSE, random_sampling parameter must be less than or equal to the smallest per-region sample size. One or more regions has a smaller number of points than the random_sampling parameter.", call. =FALSE)
      }
      region_sizes[region] <- random_sampling
    }
  }
  
  # Randomly select points from each region
  region_samples <- vector(mode = "character")
  for (region in names(region_sizes))
  {
    region_subset <- sample(rownames(stvea_object@codex_protein[stvea_object@codex_region == region,]), region_sizes[region])
    region_samples <- c(region_samples, region_subset)
  }
  # Filter dataset to only include randomly selected points from each region
  selected_points <- row.names(stvea_object@codex_protein) %in% region_samples
  if (!is.null(stvea_object@codex_protein))
  {
    stvea_object@codex_protein <- stvea_object@codex_protein[selected_points,]
  }
  if (!is.null(stvea_object@codex_size))
  {
    stvea_object@codex_size <- stvea_object@codex_size[selected_points]
  }
  if (!is.null(stvea_object@codex_cell_ID))
  {
    stvea_object@codex_cell_ID <- stvea_object@codex_cell_ID[selected_points]
  }
  if (!is.null(stvea_object@codex_region))
  {
    stvea_object@codex_region <- stvea_object@codex_region[selected_points]
  }
  if (!is.null(stvea_object@codex_blanks))
  {
    stvea_object@codex_blanks <- stvea_object@codex_blanks[selected_points,]
  }
  if (!is.null(stvea_object@codex_spatial))
  {
    stvea_object@codex_spatial <- stvea_object@codex_spatial[selected_points,]
  }
  if (!is.null(stvea_object@codex_emb))
  {
    stvea_object@codex_emb <- stvea_object@codex_emb[selected_points,]
  }
  if (!is.null(stvea_object@codex_knn))
  {
    stvea_object@codex_knn <- stvea_object@codex_knn[selected_points,]
  }
  if (!is.null(stvea_object@codex_clusters))
  {
    stvea_object@codex_clusters <- stvea_object@codex_clusters[selected_points]
  }
  if (!is.null(stvea_object@codex_clean))
  {
    stvea_object@codex_clean <- stvea_object@codex_clean[selected_points,]
  }
  
  return(stvea_object)
}

#' Create version of STvEA's FilterCODEX function for MAV-defined gating
#' (Govek et al. Science Advances 2019.) and inclusion of original region
# Function wrappers using STvEA.data class

#' Removes points from CODEX matrix that are not cells
#' as determined by the gating strategy on the blank channels
#' from the CODEX paper. Points can optionally be included/excluded by a MAV-defined vector prior to filtering.
#'
#' @param stvea_object STvEA.data class with CODEX data
#' @param size_lim lower and upper limits on size of each cell. If blank, set to 0.025 and 0.99 quantiles
#' @param blank_upper a vector with an upper bound expression cutoff for each blank channel.
#' If NULL, blank upper bounds are set as the 0.995 quantile for each blank
#' @param blank_lower a vector with a lower bound expression cutoff for each blank channel.
#' If NULL, blank lower bounds are set as the 0.002 quantile for each blank
#' @param inclusion (optional) a logical vector with an entry for each cell that indicates whether it should be included or excluded (prior to other processing).
#'
#' @export
#'
FilterCODEX <- function(stvea_object,
                              size_lim = NULL,
                              blank_upper = NULL,
                              blank_lower = NULL,
                              inclusion = NULL) {
  if (is.null(stvea_object@codex_protein) || is.null(stvea_object@codex_size) ||
      is.null(stvea_object@codex_blanks)|| is.null(stvea_object@codex_cell_ID)|| is.null(stvea_object@codex_region)) # New code added to additionally check for cell_IDs and codex_regions slots
  {
    stop("Input object must contain size, ID, and region of each CODEX cell and expression for CODEX protein and blank channels", call. =FALSE)
  }
  if (is.null(row.names(stvea_object@codex_protein))) {
    row.names(stvea_object@codex_protein) <- as.character(1:nrow(stvea_object@codex_protein))
  }
  # New code added: check that inclusion parameter is set correctly in call to function
  if (!is.null(inclusion) & (!is.logical(inclusion) | !is.vector(inclusion) | length(inclusion) != nrow(stvea_object@codex_protein)))
  {
    stop("Inclusion parameter must be a logical vector of the same length as the number of cells in the dataset.", call. =FALSE)
  }

  # New code added: filter cells by inclusion vector prior to other processing
  if(!is.null(inclusion))
  {
    stvea_object@codex_protein <- stvea_object@codex_protein[inclusion,]
    stvea_object@codex_size <- stvea_object@codex_size[inclusion]
    stvea_object@codex_blanks <- stvea_object@codex_blanks[inclusion,]
    stvea_object@codex_cell_ID <- stvea_object@codex_cell_ID[inclusion]
    stvea_object@codex_region <- stvea_object@codex_region[inclusion]
    if (!is.null(stvea_object@codex_spatial)) {
      stvea_object@codex_spatial <- stvea_object@codex_spatial[inclusion,]
    }
    if (!is.null(stvea_object@codex_clusters)) {
      stvea_object@codex_clusters <- stvea_object@codex_clusters[inclusion]
    }
    if (!is.null(stvea_object@codex_emb)) {
      stvea_object@codex_emb <- stvea_object@codex_emb[inclusion,]
    }
  }

  filter_matrix <- FilterCODEX.internal(stvea_object@codex_protein,
                                        stvea_object@codex_size,
                                        stvea_object@codex_blanks,
                                        size_lim=size_lim,
                                        blank_upper=blank_upper,
                                        blank_lower=blank_lower)
  filter <- row.names(stvea_object@codex_protein) %in% row.names(filter_matrix)
  stvea_object@codex_protein <- filter_matrix
  stvea_object@codex_size <- stvea_object@codex_size[filter]
  stvea_object@codex_blanks <- stvea_object@codex_blanks[filter,]
  stvea_object@codex_cell_ID <- stvea_object@codex_cell_ID[filter] # New code added
  stvea_object@codex_region <- stvea_object@codex_region[filter] # New code added

  if (!is.null(stvea_object@codex_spatial)) {
    stvea_object@codex_spatial <- stvea_object@codex_spatial[filter,]
  }
  if (!is.null(stvea_object@codex_clusters)) {
    stvea_object@codex_clusters <- stvea_object@codex_clusters[filter]
  }
  if (!is.null(stvea_object@codex_emb)) {
    stvea_object@codex_emb <- stvea_object@codex_emb[filter,]
  }
  return(stvea_object)
}

#' Remove region-specific artifacts caused by heterogeneity in imaging process 
#' Can remove noise-associated clusters based on size (using a fitted 
#' Gaussian distribution of size) and/or based on outlier presence in a single region 
#' @param stvea_object STvEA object
#' @param size_threshold (optional) percentile threshold for Gaussian cluster removal by size (high-pass)
#' @param outlier_threshold (optional) percentile threshold for cluster removal by region outlier (low-pass)
#'
#' @return stvea_object list of processed Monocle and STvEA objects
#'
#' @export
#'
RemoveArtifacts <- function(stvea_object,
                            size_threshold = NULL,
                            outlier_threshold = NULL) {
  # If specified, eliminate excessively small clusters below a percentile of cluster size Gaussian
  if (!is.null(size_threshold)) 
  {
    cluster_constituents <- vector(mode = "numeric", length = length(unique(stvea_object@codex_clusters)))
    names(cluster_constituents) <- as.character(unique(stvea_object@codex_clusters))
    for (cluster in names(cluster_constituents))
    {
      cluster_constituents[cluster] <- sum(stvea_object@codex_clusters == cluster)
    }
    size_threshold <- quantile(cluster_constituents,size_threshold)
    noise <- names(cluster_constituents[cluster_constituents < size_threshold])
    denoise_filter <- !(as.character(stvea_object@codex_clusters) %in% noise)
    stvea_object@codex_protein <- stvea_object@codex_protein[denoise_filter,]
    stvea_object@codex_size <- stvea_object@codex_size[denoise_filter]
    stvea_object@codex_blanks <- stvea_object@codex_blanks[denoise_filter,]
    stvea_object@codex_cell_ID <- stvea_object@codex_cell_ID[denoise_filter]
    stvea_object@codex_region <- stvea_object@codex_region[denoise_filter]
    stvea_object@codex_spatial <- stvea_object@codex_spatial[denoise_filter,]
    stvea_object@codex_clusters <- stvea_object@codex_clusters[denoise_filter]
    stvea_object@codex_clusters <- droplevels(stvea_object@codex_clusters)
    stvea_object@codex_emb <- stvea_object@codex_emb[denoise_filter,]
    stvea_object@codex_knn <- stvea_object@codex_knn[denoise_filter,]
    stvea_object@codex_clean <- stvea_object@codex_clean[denoise_filter,]
  }
  
  # If specified, remove region-specific outliers
  if (!is.null(outlier_threshold)) 
  {    
    region_clusters <- matrix(nrow = length(unique(stvea_object@codex_clusters)), ncol = length(unique(stvea_object@codex_region)))
    rownames(region_clusters) <- as.character(unique(stvea_object@codex_clusters))
    colnames(region_clusters) <- as.character(unique(stvea_object@codex_region))
    region_filters <- list()
    cluster_filters <- list()
    for (region in colnames(region_clusters))
    {
      region_filters[[region]] <- stvea_object@codex_region == region
    }
    for (cluster in rownames(region_clusters))
    {
      cluster_filters[[cluster]] <- stvea_object@codex_clusters == cluster
    }
    for (region_to_count in colnames(region_clusters))
    {
      for (cluster_to_count in rownames(region_clusters))
      {
        region_clusters[as.character(cluster_to_count),as.character(region_to_count)] <- sum(region_filters[[region_to_count]] & cluster_filters[[cluster_to_count]])
      }
    }
    region_percents <- matrix(nrow = length(unique(stvea_object@codex_clusters)), ncol = length(unique(stvea_object@codex_region)))
    rownames(region_percents) <- as.character(unique(stvea_object@codex_clusters))
    colnames(region_percents) <- as.character(unique(stvea_object@codex_region))
    for (region_to_calculate in colnames(region_percents))
    {
      for (cluster_to_calculate in rownames(region_percents))
      {
        region_percents[as.character(cluster_to_calculate),as.character(region_to_calculate)] <- region_clusters[as.character(cluster_to_calculate),as.character(region_to_calculate)]/sum(region_clusters[as.character(cluster_to_calculate),])*100
      }
    }
    outliers <- names(which(rowSums(region_percents > (outlier_threshold*100)) > 0))
    outlier_filter <- !(as.character(stvea_object@codex_clusters) %in% outliers)
    stvea_object@codex_protein <- stvea_object@codex_protein[outlier_filter,]
    stvea_object@codex_size <- stvea_object@codex_size[outlier_filter]
    stvea_object@codex_blanks <- stvea_object@codex_blanks[outlier_filter,]
    stvea_object@codex_cell_ID <- stvea_object@codex_cell_ID[outlier_filter]
    stvea_object@codex_region <- stvea_object@codex_region[outlier_filter]
    stvea_object@codex_spatial <- stvea_object@codex_spatial[outlier_filter,]
    stvea_object@codex_clusters <- stvea_object@codex_clusters[outlier_filter]
    stvea_object@codex_clusters <- droplevels(stvea_object@codex_clusters)
    stvea_object@codex_emb <- stvea_object@codex_emb[outlier_filter,]
    stvea_object@codex_knn <- stvea_object@codex_knn[outlier_filter,]
    stvea_object@codex_clean <- stvea_object@codex_clean[outlier_filter,]
  }
  
  return(stvea_object)
}

#' Create version of STvEA's PlotExprCODEXspatial function that allows for
#' plotting of UMAP clusters and specific regions
#' Plot expression of a gene or protein in the CODEX spatial coordinates
#'
#' @param stvea_object STvEA.data class with CODEX expression and spatial xy
#' @param region name of region to plot
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA", "protein", or "cluster"
#' @param high_color color of high expression on the color ramp
#' @param high_color2 color of high expression for second protein or gene, if provided
#' @param low_color color of low/no expression on the color ramp
#' @param cluster_colors vector indicating color of each cluster, in the format c(cluster1 = "blue", cluster2 = "red", etc.)
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotExprCODEXspatial <- function(stvea_object,
                                       region = NULL,
                                       name,
                                       type="protein",
                                       high_color="red",
                                       high_color2="green",
                                       low_color="white",
                                       cluster_colors=NULL,
                                       pt_size =0.8) {
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    } else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    } else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
  } else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  } else if (type == "cluster") {
    plotting_data <- as.character(stvea_object@codex_clusters) # New code added: plotting of clusters
  } else {
    stop("type must be \"RNA\", \"protein\", or \"cluster\".", call. =FALSE)
  }

  # New code added: restrict plotting data to region specified by user
  if (!is.null(region) & (!is.character(region) | !is.vector(region) | length(region) != 1))
  {
    stop("Region parameter must be a single character.", call. =FALSE)
  }
  if (!is.null(region))
  {
    if (type != "cluster")
    {
      plotting_data <- plotting_data[as.character(stvea_object@codex_region) == region,]
    } else
    {
      plotting_data <- plotting_data[as.character(stvea_object@codex_region) == region]
    }
    stvea_object@codex_spatial <- stvea_object@codex_spatial[as.character(stvea_object@codex_region) == region,]
  }

  # New code added: do not set a continuous color palette for cluster plotting
  if (type != "cluster")
  {
    if (length(name) > 2) {
      stop("name must be at most length 2", call. =FALSE)
    }

    rbPal1 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color,1)), alpha=TRUE)
    color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]
    subtitle <- paste("Expression of ", name[1], " (", high_color, ")", sep="")

    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color,0),alpha(high_color2,1)), alpha=TRUE)
      color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep="")
    }
  } else
  {
    subtitle <- "Cluster visualization"
  }


  x_tmp <- stvea_object@codex_spatial[,"x"]
  x_tmp <- x_tmp - min(x_tmp)
  y_tmp <- stvea_object@codex_spatial[,"y"]
  y_tmp <- y_tmp - min(y_tmp)
  spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))

  # New code added: generate distinct plot with legend for cluster plotting
  if (type != "cluster")
  {
    ggplot(spatial_tmp,
           aes(x=x,y=y,color=factor(1:length(color)))) +
      geom_point(size=pt_size, alpha=0.5) +
      scale_color_manual(values=alpha(color,1)) +
      guides(color=FALSE) +
      ylim(max(y_tmp), 0) +
      labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
      theme_void() + coord_fixed()
  } else if (is.null(cluster_colors))
  {
    ggplot(spatial_tmp,
           aes(x=x,y=y,color=plotting_data)) +
      geom_point(size=pt_size, alpha=0.5) +
      guides(colour=guide_legend(title='UMAP Cluster', override.aes = list(size=5))) +
      ylim(max(y_tmp), 0) +
      labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
      theme_void() + coord_fixed()
  } else
  {
    ggplot(spatial_tmp,
           aes(x=x,y=y,color=plotting_data)) +
      geom_point(size=pt_size, alpha=0.5) +
      guides(colour=guide_legend(title='UMAP Cluster', override.aes = list(size=5))) +
      ylim(max(y_tmp), 0) +
      labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
      theme_void() + coord_fixed() + scale_colour_manual(values = cluster_colors)
  }
}

#' Given a list of coordinates for n datapoints, identify the neighbors of each datapoint
#' within a specified search distance and maximum number of neighbors.
#' Returns an n x n adjacency matrix of 0's/1's to indicate absent/present neighbor relationships.
#'
#' @param coords a n x m matrix of n datapoints with m coordinates
#' @param n_neighbors maximum number of neighbors to identify per datapoint (by default, all within range)
#' @param dist distance to search for neighbors for each datapoint
#'
#' @export
#'

dist_graph <- function(coords, n_neighbors, dist) {
  n_points <- nrow(coords)
  # To prevent integer overflow in RANN, subsample excessive size requests to maximum
  if (n_points*as.numeric(n_neighbors) >= (2^31))
  {
    coords <- coords[sample(1:n_points, size = ceiling(2^31/n_neighbors-1), replace=F),]
    if (n_neighbors > nrow(coords))
      n_neighbors <- nrow(coords)
  }
  # Get indices of neighboring cells within specified pixel distance
  nn_idx <- RANN::nn2(coords, searchtype = "radius", k=n_neighbors, radius=dist)$nn.idx 
  # Convert matrix of neighbor indices to a two-column list/matrix of neighbor-neighbor connections
  neighbor_connections <- cbind(nn_idx[, 1], c(nn_idx[, -1]))
  # Remove null neighbor-neighbor connections
  neighbor_connections <- neighbor_connections[neighbor_connections[,2] != 0,]
  # If no neighbor-neighbor connections are left, return NULL
  if (nrow(neighbor_connections) == 0)
  {
    return(NULL)
  } else
  {
    # If the last indexed cell doesn't have a neighbor-neighbor connection, populate it with a dummy connection to prevent igraph from truncating adjacency matrix
    if (max(neighbor_connections) < n_points)
    {
      neighbor_connections <- rbind(neighbor_connections,c(n_points,n_points))
    }
    # Use igraph to convert neighbor-neighbor connection list to adjacency matrix
    neighbor_connections <- igraph::graph.edgelist(neighbor_connections)
    adjacency_matrix <- igraph::get.adjacency(neighbor_connections)
    adjacency_matrix <- 1*((adjacency_matrix + t(adjacency_matrix)) > 0)
    return(adjacency_matrix)
  }
}

#' For each cell phenotype pairing of interest (1 and 2), 
#' calculates either a summed interaction score or individual scores for each datapoint
#' using the dot product of feature vector f1, adjacency matrix j, and feature vector f2.
#' In this case, f1 and f2 are feature vectors of cell phenotype probability in each cell/datapoint,
#' while j is an adjacency matrix identifying neighboring datapoints within k nearest neighbors
#' (must be generated a priori using a knn graph function).
#'
#' @param adjacency_matrix adjacency matrix/graph of k-nearest-neighbors for each datapoint (n datapoints x n datapoints)
#' @param cell_matrix feature matrix of probabilities of phenotypic representation in each datapoint (m phenotypes x n datapoints, rownames must be phenotypes)
#' @param cell_pairs matrix of all phenotype pairings of interest (m phenotypes x 2 columns)
#' @param cell_labels optional vector of cell ID labels for each datapoint (n datapoints)
#' @param calculate_individual_scores boolean to calculate individual datapoint-level scores (default: FALSE)
#' @param normalize boolean to normalize scores by dataset size, i.e. number of datapoints (default: FALSE)
#'
#' @export
#' 

CalculateInteractionScores <- function(adjacency_matrix, cell_matrix, cell_pairs, cell_labels = NULL, calculate_individual_scores = FALSE, normalize = FALSE) {
  
  # If cell labels not provided, just use sequential numbers
  if (is.null(cell_labels))
  {
    cell_labels <- paste("Cell",1:nrow(adjacency_matrix), sep = "_")
  }
  
  # Generate summary data frame for later export 
  score_summary <- data.frame(f1 = cell_pairs[,1], f2 = cell_pairs[,2])
  rownames(score_summary) <- paste(cell_pairs[,1], cell_pairs[,2], sep = " - ")
  
  # Calculate datapoint-level interaction score for all possible cell pairs
  if (calculate_individual_scores)
  {
    sample_scores <- matrix(nrow = nrow(cell_pairs), ncol = length(cell_labels))
    colnames(sample_scores) <- cell_labels
    for (i in 1:nrow(cell_pairs)) 
    {
      f1 <- cell_matrix[cell_pairs[i,1],]
      f2 <- cell_matrix[cell_pairs[i,2],]
      sample_scores[i,] <- as.vector(f1%*%adjacency_matrix*f2)
    }
    
    # Populate scores and export data frame
    for (datapoint in cell_labels)
    {
      score_summary[[datapoint]] <- sample_scores[,datapoint]
    }
  }
  
  # Calculate sample-level interaction score for all possible cell pairs
  if (!calculate_individual_scores)
  {
    sample_scores <- vector(mode = "numeric", length = nrow(cell_pairs))
    for (i in 1:nrow(cell_pairs)) 
    {
      f1 <- cell_matrix[cell_pairs[i,1],]
      f2 <- cell_matrix[cell_pairs[i,2],]
      sample_scores[i] <- sum(f1%*%adjacency_matrix*f2)
    }
    
    # Populate scores and export data frame
    score_summary$score <- sample_scores
  }
  
  # If specified, normalize scores by dataset size (number of datapoints)
  if (normalize)
  {
    for (score_category in colnames(score_summary)[!colnames(score_summary) %in% c("f1","f2")])
    {
      score_summary[[score_category]] <- score_summary[[score_category]]/nrow(adjacency_matrix)*1000
    }
  }
  
  # Export data frame
  return(score_summary)
}


#' Create version of STvEA's heatmap function that uses ggplot2's heatmap instead
#' Visualize heatmap of Adjacency Score run on
#' all combinations of feature pairs
#'
#' @param adj_score_output output of any of the AdjScore functions
#' - a matrix where each row contains the features in a pair and their
#' Adjacency Score and q values
#' @param score_parameter score parameter to plot in heatmap (either "score" (default), "p", "q", or "log10q")
#'
#' @export
#'
AdjScoreHeatmap <- function(adj_score_output, score_parameter = "score") {
  # New code added: check that score_parameter is set correctly
  if (!is.null(score_parameter) & (score_parameter != "score" & score_parameter != "p" & score_parameter != "q" & score_parameter != "log10q"))
  {
    stop("score_parameter must be either \"score\", \"p\", \"q\", or \"log10q\".", call. =FALSE)
  }
  # New code added: create appropriately structured dataframe for heatmap plotting
  adj_score_output <- cbind(adj_score_output,matrix(ncol = 1, nrow = nrow(adj_score_output)))
  colnames(adj_score_output)[6] <- "log10q"
  adj_score_output[,"log10q"] <- log10(adj_score_output[,"q"])
  temp <- adj_score_output
  temp[,"f1"] <- adj_score_output[,"f2"]
  temp[,"f2"] <- adj_score_output[,"f1"]
  adj_score_output <- rbind(adj_score_output,temp)
  score_data_frame <- as.data.frame(adj_score_output)
  # Set NA values to minimum for plotting
  score_data_frame[is.infinite(score_data_frame[,"log10q"]),"log10q"] <- min(score_data_frame[is.finite(score_data_frame[,"log10q"]),"log10q"])
  # Set legend titles according to score parameter used
  if (score_parameter == "score")
  {
    legend_title <- "Raw Score"
  }
  if (score_parameter == "p")
  {
    legend_title <- "p Value"
  }
  if (score_parameter == "q")
  {
    legend_title <- "Benjamini-Hochberg q Value"
  }
  if (score_parameter == "log10q")
  {
    legend_title <- "Log10(q) Value"
  }
  # Plot heatmap using ggplot2's heatmap functionality
  ggplot(score_data_frame, aes(f, g, fill = eval(as.name(score_parameter)))) + geom_tile() + scale_fill_viridis_c() + guides(fill=guide_colourbar(title = legend_title, ticks = FALSE))
}

#' Create version of STvEA's heatmap function that outputs adjacency matrix
#' Output matrix of Adjacency Score run on
#' all combinations of feature pairs
#'
#' @param adj_score_output output of any of the AdjScore functions
#' - a matrix where each row contains the features in a pair and their
#' Adjacency Score and q values
#' @param score_parameter score parameter to output in adjacency matrix (either "score" (default), "p", "q", or "log10q")
#' @param top_interactions restricts the output to the top interactions in the dataset. Must be a number indicating the number of interactions to keep. Does not keep self-interactions.
#' @param decreasing if true (default), sorts top interactions in decreasing order (i.e. highest numbers are kept).
#'
#' @export
#'
AdjScoreHeatmapMatrix <- function(adj_score_output, score_parameter = "score", top_interactions = NULL, decreasing = TRUE) {
  # New code added: check that score_parameter and top_interactions are set correctly
  if (!is.null(score_parameter) & (score_parameter != "score" & score_parameter != "p" & score_parameter != "q" & score_parameter != "log10q"))
  {
    stop("score_parameter must be either \"score\", \"p\", \"q\", or \"log10q\".", call. =FALSE)
  }
  if (!is.null(top_interactions) & (!is.numeric(top_interactions) | !is.vector(top_interactions) | length(top_interactions) != 1))
  {
    stop("Top_interactions parameter must be a single number.", call. =FALSE)
  }
  if (!is.null(top_interactions))
  {
    top_interactions <- round(top_interactions)
  }

  if (!is.null(top_interactions) & decreasing)
  {
    adj_score_output <- adj_score_output[order(abs(as.numeric(adj_score_output[,score_parameter])), decreasing=TRUE),]
  } else if (!is.null(top_interactions))
  {
    adj_score_output <- adj_score_output[order(abs(as.numeric(adj_score_output[,score_parameter])), decreasing=FALSE),]
  }

  if (!is.null(top_interactions))
  {
    adj_score_output <- adj_score_output[adj_score_output[,"f1"] != adj_score_output[,"f2"],]
    adj_score_output <- adj_score_output[1:min(top_interactions,nrow(adj_score_output)),]
  }

  nodes <- unique(c(unique(adj_score_output[,"f1"]),unique(adj_score_output[,"f2"])))
  nodes <- sort(nodes)
  heatmap_matrix <- matrix(0, ncol = length(nodes), nrow = length(nodes), dimnames = list(nodes,nodes))

  if (score_parameter == "log10q")
  {
    for (i in 1:nrow(adj_score_output))
    {
      heatmap_matrix[as.character(adj_score_output[i,"f1"]),as.character(adj_score_output[i,"f2"])] <- log10(adj_score_output[i,"q"]+1e-15)
      heatmap_matrix[as.character(adj_score_output[i,"f2"]),as.character(adj_score_output[i,"f1"])] <- log10(adj_score_output[i,"q"]+1e-15)
    }
  } else
  {
    for (i in 1:nrow(adj_score_output))
    {
      heatmap_matrix[as.character(adj_score_output[i,"f1"]),as.character(adj_score_output[i,"f2"])] <- adj_score_output[i,score_parameter]
      heatmap_matrix[as.character(adj_score_output[i,"f2"]),as.character(adj_score_output[i,"f1"])] <- adj_score_output[i,score_parameter]
    }
  }
  return(heatmap_matrix)
}

#' Create version of STvEA's PlotExprCODEXspatial function that allows for
#' plotting of UMAP clusters and specific regions
#' Plot expression of a gene or protein in the CODEX spatial coordinates
#'
#' @param stvea_object STvEA.data class with CODEX expression and spatial xy
#' @param region name of region to plot
#' @param name gene or protein name to plot
#' @param type which type of expression data should be plotted? "RNA", "protein", or "cluster"
#' @param high_color color of high expression on the color ramp
#' @param high_color2 color of high expression for second protein or gene, if provided
#' @param low_color color of low/no expression on the color ramp
#' @param cluster_colors vector indicating color of each cluster, in the format c(cluster1 = "blue", cluster2 = "red", etc.)
#' @param pt_size ggplot geom_point size
#'
#' @export
#'
PlotExprCODEXspatialNoAlpha <- function(stvea_object,
                                        region = NULL,
                                        name,
                                        type="protein",
                                        high_color="red",
                                        high_color2="green",
                                        low_color="white",
                                        cluster_colors=NULL,
                                        pt_size =0.8) {
  if (is.null(stvea_object@codex_spatial)) {
    stop("stvea_object does not contain CODEX spatial information")
  }
  if (type == "protein") {
    if (!is.null(stvea_object@codex_clean)) {
      plotting_data <- stvea_object@codex_clean
    } else if (!is.null(stvea_object@codex_protein)) {
      plotting_data <- stvea_object@codex_protein
    } else {
      stop("stvea_object must contain CODEX protein data with type=\"protein\"")
    }
  } else if (type == "RNA") {
    plotting_data <- stvea_object@codex_mRNA
  } else if (type == "cluster") {
    plotting_data <- as.character(stvea_object@codex_clusters) # New code added: plotting of clusters
  } else {
    stop("type must be \"RNA\", \"protein\", or \"cluster\".", call. =FALSE)
  }
  
  # New code added: restrict plotting data to region specified by user
  if (!is.null(region) & (!is.character(region) | !is.vector(region) | length(region) != 1))
  {
    stop("Region parameter must be a single character.", call. =FALSE)
  }
  if (!is.null(region))
  {
    if (type != "cluster")
    {
      plotting_data <- plotting_data[as.character(stvea_object@codex_region) == region,]
    } else
    {
      plotting_data <- plotting_data[as.character(stvea_object@codex_region) == region]
    }
    stvea_object@codex_spatial <- stvea_object@codex_spatial[as.character(stvea_object@codex_region) == region,]
  }
  
  # New code added: do not set a continuous color palette for cluster plotting
  if (type != "cluster")
  {
    if (length(name) > 2) {
      stop("name must be at most length 2", call. =FALSE)
    }
    
    rbPal1 <- colorRampPalette(c(alpha(low_color,1),alpha(high_color,1)), alpha=TRUE)
    color <- rbPal1(100)[as.numeric(cut(plotting_data[,name[1]],breaks = 100))]
    subtitle <- paste("Expression of ", name[1], " (", high_color, ")", sep="")
    
    if (length(name) == 2) {
      rbPal2 <- colorRampPalette(c(alpha(low_color,1),alpha(high_color2,1)), alpha=TRUE)
      color2 <- rbPal2(100)[as.numeric(cut(plotting_data[,name[2]],breaks = 100))]
      color <- sapply(1:length(color), function(m) colorRampPalette(c(color[m],color2[m]), alpha=TRUE)(3)[2])
      subtitle <- paste(subtitle, " and ", name[2], " (", high_color2, ")", sep="")
    }
  } else
  {
    subtitle <- "Cluster visualization"
  }
  
  
  x_tmp <- stvea_object@codex_spatial[,"x"]
  x_tmp <- x_tmp - min(x_tmp)
  y_tmp <- stvea_object@codex_spatial[,"y"]
  y_tmp <- y_tmp - min(y_tmp)
  spatial_tmp <- as.data.frame(cbind(x = x_tmp, y = y_tmp))
  
  # New code added: generate distinct plot with legend for cluster plotting
  if (type != "cluster")
  {
    ggplot(spatial_tmp,
           aes(x=x,y=y,color=factor(1:length(color)))) +
      geom_point(size=pt_size) +
      scale_color_manual(values=alpha(color,1)) +
      guides(color=FALSE) +
      ylim(max(y_tmp), 0) +
      labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
      theme_void() + coord_fixed()
  } else if (is.null(cluster_colors))
  {
    ggplot(spatial_tmp,
           aes(x=x,y=y,color=plotting_data)) +
      geom_point(size=pt_size) +
      guides(colour=guide_legend(title='UMAP Cluster', override.aes = list(size=5))) +
      ylim(max(y_tmp), 0) +
      labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
      theme_void() + coord_fixed()
  } else
  {
    ggplot(spatial_tmp,
           aes(x=x,y=y,color=plotting_data)) +
      geom_point(size=pt_size) +
      guides(colour=guide_legend(title='UMAP Cluster', override.aes = list(size=5))) +
      ylim(max(y_tmp), 0) +
      labs(title = paste("Spatial",type,"expression"), subtitle = subtitle) +
      theme_void() + coord_fixed() + scale_colour_manual(values = cluster_colors)
  }
}