adjacencyMatrix <- function(clones,
                            dist_type = "levenshtein",
                            max_dist = 1,
                            sparse = TRUE) {

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
  if (dist_type == "levenshtein") {
    if (sparse) {
      out <- levAdjacencyMatSparse(clones, max_dist)
    } else {
      out <- levAdjacencyMatDense(clones, max_dist) }
  }

  return(out)

}