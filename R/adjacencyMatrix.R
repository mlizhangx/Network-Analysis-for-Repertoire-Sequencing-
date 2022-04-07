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
  if (dist_type == "levenshtein") {
      out <- levAdjacencyMatSparse(clones, max_dist)
  } else {
    stop('"levenshtein" is currently the only available option for `dist_type`')
  }


  return(out)

}