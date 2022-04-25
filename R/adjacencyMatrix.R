adjacencyMatrix <- function(clones,
                            dist_type = "levenshtein",
                            max_dist = 1) {

  # attempt to coerce clones to character vector
  if (length(clones) == 0) stop("'clones' has zero length")
  clones <- as.vector(clones, mode = "character")
  if (!is.character(clones)) {
    stop("'clones' must be cocercible to a character vector")
  }
  if (!is.vector(clones)) {
    stop("'clones' must be cocercible to a character vector")
  }

  # Compute adjacency matrix
  if (dist_type %in% c("levenshtein", "Levenshtein, lev, Lev, l, L")) {
    out <- levAdjacencyMatSparse(clones, max_dist)
  } else if (dist_type %in% c("hamming", "Hamming", "ham", "Ham", "h", "H")) {
    out <- hamAdjacencyMatSparse(clones, max_dist)
  } else {
    stop('supported options for `dist_type` include "levenshtein" and "hamming"')
  }

  # Number of nodes with positive network degree
  num_nodes <- dim(out)[[1]]

  if (num_nodes == 0) {

    warning(paste0(
      "Network contains no edges; no pairs of sequences are closer than the specified distance threshold (at most ",
      max_dist, ") for adjacency. Try using a higher threshold using the `max_dist` argument"))

  } else {

    cat(paste0(
      "Adjacency matrix created for the ", num_nodes,
      " network nodes with positive degree, i.e., all sequences that are within the specified distance threshold (at most ",
      max_dist, ") of at least one other sequence.\n"))
    # Import record of selected column IDs and use for matrix row names
    clone_ids <- utils::read.table("col_ids.txt")
    dimnames(out)[[1]] <- clone_ids$V1
    dimnames(out)[[2]] <- clones[clone_ids$V1]
    cat(paste0(
      "The row names of the adjacency matrix contain the original index values of the corresponding sequences; the column names contain the sequences themselves. They can be accessed using `dimnames()`\n"))
  }

  # Remove temporary file of column ids
  file.remove("col_ids.txt")

  return(out)

}
