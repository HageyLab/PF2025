## Previously developed in Mascharak & Guo et al. Cell Reports Medicine 2023.

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
#' @param log_scale boolean to log-scale scores (default: FALSE)
#'
#' @export
#' 

CalculateInteractionScores <- function(adjacency_matrix, cell_matrix, cell_pairs, cell_labels = NULL, calculate_individual_scores = FALSE, normalize = FALSE, log_scale = FALSE) {
  
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
  
  # If specified, take log of scores
  if (log_scale)
  {
    for (score_category in colnames(score_summary)[!colnames(score_summary) %in% c("f1","f2")])
    {
      score_summary[[score_category]][score_summary[[score_category]] != 0] <- log(score_summary[[score_category]][score_summary[[score_category]] != 0]*1000)
    }
  }
  
  # Export data frame
  return(score_summary)
}


#' Create version of STvEA's heatmap function that outputs adjacency matrix (Govek et al. Science Advances 2019.)
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


#' Basic function to convert individual rows of a matrix to z-scores
#'
#' @param mat matrix to convert
#' @param center center values around mean (default: TRUE)
#' @param scale scale values by standard deviation (default: TRUE)
#'
#' @export
#'
scale_rows <- function(mat, center = TRUE, scale = TRUE) {
  # Calculate means at row-level, if specified
  if (center)
  {
    means <- rowMeans(mat, na.rm = TRUE)
  } else
  {
    means <- rep(0, length = nrow(mat))
  }
  # Calculate standard deviations at row-level, if specified
  if (scale) 
  {
    stdevs <- matrixStats::rowSds(mat, center = means)
  } else {
    stdevs <- rep(1, length = nrow(mat))
  }
  mat <- (mat - means) / stdevs
  return(mat)
}
