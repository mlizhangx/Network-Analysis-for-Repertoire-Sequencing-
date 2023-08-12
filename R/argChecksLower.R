
# Universal Checks ------------------------------------------------------------

.stopifnot <- function(condition, argname, message) {
  if (!condition) {
    stop(paste(argname, message))
  }
}

.nonNull <- function(value, argname) {
  if (is.null(value)) {
    stop(paste(argname, "is required but value is NULL"))
  }
}

.noNAs <- function(value, argname) {
  if (sum(is.na(value)) > 0) {
    stop(paste(argname, "contains NAs"))
  }
}

.hasLength1 <- function(value, argname) {
  if (length(value) != 1) {
    stop(paste(argname, "must have length 1"))
  }
}

.hasLength2 <- function(value, argname) {
  if (length(value) != 2) {
    stop(paste(argname, "must have length 2"))
  }
}

.hasLength <- function(length_spec, value, argname) {
  if (length(value) != length_spec) {
    stop(paste(argname, "must have length", length_spec))
  }
}

.hasPosLength <- function(value, argname) {
  if (length(value) == 0) {
    stop(paste(argname, "must have positive length"))
  }
}



# Type Checks -------------------------------------------------------------

# checks that an argument is strictly TRUE or FALSE
.isTF <- function(value, argname) {
  .hasLength1(value, argname)
  .stopifnot(sum(value %in% c(TRUE, FALSE)) == 1,
             argname, "must evaluate to \"TRUE\" or \"FALSE\""
  )
}

.isTFOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isTF(value, argname)
  }
}

.isTFOrAuto <- function(value, argname) {
  .nonNull(value, argname)
  .hasLength1(value, argname)
  if (is.character(value)) {
    .stopifnot(
      value == "auto",
      argname, "must be TRUE, FALSE, or 'auto'1"
    )
  } else {
    .stopifnot(
      sum(value %in% c(TRUE, FALSE)) == 1,
      argname, "must be TRUE, FALSE, or 'auto'2"
    )
  }
}

.isNumeric <- function(value, argname) {
  .stopifnot(
    is.numeric(value),
    argname, "must be of type 'numeric'"
  )
}

.isNumericOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isNumeric(value, argname)
  }
}

.isFinite <- function(value, argname) {
  .stopifnot(
    all(is.finite(value)),
    argname, "contains non-finite values"
  )
}

# checks that argument is nonnegative and finite
.isNonneg0 <- function(value, argname) {
  .hasLength1(value, argname)
  .isNumeric(value, argname)
  .isFinite(value, argname)
  .stopifnot(
    value >= 0,
    argname, "must be nonnegative"
  )
}

.isNonneg <- function(value, argname) {
  .nonNull(value, argname)
  .isNonneg0(value, argname)
}
.isNonnegOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isNonneg0(value, argname)
  }
}

# checks that argument is strictly positive and finite
.isPos0 <- function(value, argname) {
  .hasLength1(value, argname)
  .isNumeric(value, argname)
  .isFinite(value, argname)
  .stopifnot(
    value > 0,
    argname, "must be strictly positive"
  )
}
.isPos <- function(value, argname) {
  .nonNull(value, argname)
  .isPos0(value, argname)
}
.isPosOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isPos0(value, argname)
  }
}

.isInt <- function(value, argname) {
  .stopifnot(value %% 1 == 0,
             argname, "must be an integer"
  )
}

# checks that argument is a strictly positive integer
.isPosInt0 <- function(value, argname) {
  .isPos(value, argname)
  .isInt(value, argname)
}

.isPosInt <- function(value, argname) {
  .nonNull(value, argname)
  .isPosInt0(value, argname)
}
.isPosIntOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isPosInt0(value, argname)
  }
}

# check of type character
.isChar <- function(value, argname) {
  .stopifnot(
    is.character(value), argname,
    "must be of type 'character'"
  )
}

# is a string
.isString0 <- function(value, argname) {
  .isChar(value, argname)
  .hasLength1(value, argname)
}

# is a string and is not NULL
.isString <- function(value, argname) {
  .nonNull(value, argname)
  .isString0(value, argname)
}

# is a string or is NULL
.isStringOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isString0(value, argname)
  }
}


# check of type character or numeric
.isCharOrNumeric <- function(value, argname) {
  .stopifnot(
    is.character(value) || is.numeric(value),
    argname, "must be of type 'character' or 'numeric'"
  )
}


# check of type character or numeric
.isCharOrNumericOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isCharOrNumeric(value, argname)
  }
}

.isCharOrNumericScalar <- function(value, argname) {
  .isCharOrNumeric(value, argname)
  .hasLength1(value, argname)
  if (is.numeric(value)) {
    .isFinite(value, argname)
  }
}

.isCharOrNumericScalarOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isCharOrNumericScalar(value, argname)
  }
}

.isCharVector <- function(value, argname) {
  .stopifnot(
    is.vector(value, mode = "character"),
    argname, "must be a character vector"
  )
  .hasPosLength(value, argname)
}

.isCharVectorOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isCharVector(value, argname)
  }
}

.isCharOrNumericVector <- function(value, argname) {
  .stopifnot(
    is.vector(value, mode = "character") || is.vector(value, mode = "numeric"),
    argname, "must be a character or numeric vector"
  )
  if (is.numeric(value)) {
    .isFinite(value, argname)
  }
}

.isCharOrNumericVectorOrNull <- function(value, argname) {
  if (!is.null(value)) {
    .isCharOrNumericVector(value, argname)
  }
}

.isIgraph <- function(value, argname) {
  .stopifnot(
    inherits(value, "igraph"),
    argname, "must be of class igraph"
  )
}

.isGgraph <- function(value, argname) {
  .stopifnot(
    inherits(value, "ggraph"),
    argname, "must be of class ggraph"
  )
}

.isDataFrame <- function(value, argname) {
  .stopifnot(
    is.data.frame(value),
    argname, "must be a data frame"
  )
}

.isList <- function(value, argname) {
  .stopifnot(
    is.list(value),
    argname, "must be a list"
  )
}

.isAdjacencyMatrix <- function(value, argname) {
  .stopifnot(
    is.matrix(value) || inherits(value, "dgCMatrix"),
    argname, "must be a matrix or dgCMatrix"
  )
  .stopifnot(
    dim(value)[[1]] == dim(value)[[2]],
    argname, "must have the same row and column dimensions"
  )
  .stopifnot(
    all(unique(as.vector(value)) %in% c(0, 1)),
    argname, "has entries with values other than 0 or 1"
  )
}

# Data Frames and Lists ----------------------------------------------------


.hasAtLeastTwoRows <- function(data) {
  # has at least two rows
  stopifnot("need at least two data rows" = nrow(data) > 1)
}


# check that argument is a valid reference for a single column of `data`
.isDataCol0 <- function(data, colref, argname) {
  .hasLength1(colref, argname)
  if (is.character(colref)) {
    .stopifnot(
      colref %in% names(data),
      argname, "does not specify a valid column name within the data")
  } else if (is.numeric(colref)) {
    .stopifnot(
      colref >= 1 && colref <= ncol(data),
      argname, "specifies a column number outside of the valid range")
  } else {
    stop(
      paste(argname, "value has type other than 'character' or 'numeric'")
    )
  }
}

.isDataCol <- function(data, colref, argname) {
  .nonNull(colref, argname)
  .isDataCol0(data, colref, argname)
}

# check that argument is a valid reference for a single column of `data` or NULL
.isDataColOrNull <- function(data, colref, argname) {
  if (!is.null(colref)) { .isDataCol0(data, colref, argname) }
}

# check that argument is a valid reference for one or more columns of `data`, or NULL
.isDataColsOrNull <- function(data, colref, argname) {
  if (!is.null(colref)) {
    .hasPosLength(colref, argname)
    if (length(colref) == 1) {
      .isDataCol0(data, colref, argname)
    } else if (length(colref) > 1) {
      for (i in length(colref)) {
        .isDataCol0(data, colref[[i]], argname)
      }
    }
  }
}

# Checks that column specified for sequence column is valid
.isSeqCol <- function(data, colref) {
  .nonNull(colref, "`seq_col`")
  if (length(colref) == 1) {
    .isDataCol0(data, colref, "`seq_col`")
    .isValidSeqVector(data[[colref]])
  } else if (length(colref) == 2) {
    .isDataCol0(data, colref[[1]], "first entry of `seq_col`")
    .isDataCol0(data, colref[[2]], "second entry of `seq_col`")
    .isValidSeqVector(data[[colref[[1]]]])
    .isValidSeqVector(data[[colref[[2]]]])
  } else {
    stop("`seq_col` must have length 1 or 2")
  }
}

.isValidSeqVector <- function(value) {
  value <- as.vector(value, mode = "character")
  if (!is.vector(value, mode = "character")) {
    stop("specified column or vector of receptor sequences must be coercible to a character vector")
  }
}

.hasListElement <- function(value, argname, element_name) {
  .stopifnot(
    element_name %in% names(value),
    argname, paste("does not contain an element named", element_name)
  )
}

# Network Output ----------------------------------------------------------

.hasNodeAndClusterData <- function(value, argname) {
  .isList(value, argname)
  .hasPosLength(value, argname)
  .hasListElement(value, argname, "node_data")
  .hasListElement(value, argname, "cluster_data")
  .isDataFrame(value$node_data, paste0(argname, "$node_data"))
  .isDataFrame(value$cluster_data, paste0(argname, "$cluster_data"))
}

.isBaseNetworkOutput <- function(value, argname) {
  .isList(value, argname)
  .hasPosLength(value, argname)
  .hasListElement(value, argname, "node_data")
  .hasListElement(value, argname, "igraph")
  .hasListElement(value, argname, "adjacency_matrix")
  .isDataFrame(value$node_data, paste0(argname, "$node_data"))
  .isIgraph(value$igraph, paste0(argname, "$igraph"))
  .isAdjacencyMatrix(
    value$adjacency_matrix, paste0(argname, "$adjacency_matrix")
  )
  .checkIgraphAgainstData(value$igraph, value$node_data)
  .checkIgraphAgainstMatrix(value$igraph, value$adjacency_matrix)
}


.isPlotlist <- function(value, argname) {
  .isList(value, argname)
  .hasPosLength(value, argname)
  for (i in 1:length(value)) {
    .isGgraph(value[[i]], paste("each element of", argname))
  }
}

# Specific Arguments -----------------------------------------------------------




.isDistType <- function(value) {
  .isString(value, "dist_type")
  stopifnot(
    "Invalid option for `dist_type` argument" =
      value %in% c(
        "levenshtein", "Levenshtein, lev, Lev, l, L",
        "hamming", "Hamming", "ham", "Ham", "h", "H"
      )
  )
}



.isInputType <- function(value) {
  .isString(value, "input_type")
  valid_input_types <- c("csv", "table", "tsv", "txt", "rds", "rda")
  if (!value %in% valid_input_types) {
    valid_input_types <- paste(valid_input_types, collapse = ", ")
    stop(
      paste("input_type must be one of:", valid_input_types)
    )
  }
}

# type = "network" for SaveNetwork()
# type = "generic" for .saveDataGeneric(), used in findAssociatedClones()
.isOutputType <- function(value, type = "network") {
  .isString(value, "output_type")
  valid_types <- switch(
    type,
    "network" = c("individual", "rds", "rda"),
    "findPublicClusters" = c("rds", "rda", "csv"),
    "findAssociatedClones" = c("csv", "tsv", "rds", "rda"),
    "generic" = c("rds", "rda", "csv", "tsv", "table")
  )
  if (!value %in% valid_types) {
    default_type <- switch(
      type,
      "network" = "rda",
      "findPublicClusters" = "rds",
      "findAssociatedClones" = "csv",
      "generic" = "rda"
    )
    warning(
      paste("output_type is invalid. Defaulting to", default_type)
    )
  }
}

# checks that value is NULL, "auto", or valid column reference
.checkColorNodesBySingle <- function(value, data, node_stats) {
  valid_names <- names(data)
  if (node_stats) {
    valid_names <- c(valid_names, names(chooseNodeStats()))
    valid_names <- valid_names[-length(valid_names)]  # remove "all_stats"
  }
  stopifnot(
    "`color_nodes_by` specifies one or more variables not present in `data` or among the node-level network properties to be computed" =
      value %in% valid_names
  )
}

.checkColorNodesBy <- function(value, data, node_stats = FALSE, plots = TRUE) {
  if (!is.null(value) && isTRUE(plots)) {
    .isCharVector(value, "color_nodes_by")
    if (length(value) == 1) {
      if (value != "auto") {
        .checkColorNodesBySingle(value, data, node_stats)
      }
    } else { # multiple values; check each individually
      for (i in seq_along(value)) {
        .checkColorNodesBySingle(value[[i]], data, node_stats)
      }
    }
  }
}


.checkColorSchemeSingle <- function(value) {
  valid_names <- c(
    "default", "viridis", "magma", "inferno", "plasma", "cividis", "rocket",
    "mako", "turbo", "A", "B", "C", "D", "E", "F", "G", "H"
  )
  valid_names <- c(
    valid_names,
    paste0(valid_names, "-1"),
    grDevices::hcl.pals()
  )
  stopifnot(
    "`color_scheme` contains one or more values which are not supported (see ?plotNetworkGraph)" =
      value %in% valid_names
  )
}


.checkColorScheme <- function(value, color_nodes_by, plots = TRUE) {
  if (isTRUE(plots) && !is.null(color_nodes_by) && !is.null(value)) {
    .isCharVector(value, "color_scheme")
    .stopifnot(
      length(value) == 1 ||
        length(value) == length(color_nodes_by),
      "color_scheme", "must have length 1 or the same length as color_nodes_by"
    )
    for (i in seq_along(value)) {
      .checkColorSchemeSingle(value[[i]])
    }
  }

}

.checkColorTitle <- function(value, color_nodes_by) {
  if (!is.null(color_nodes_by)) {
    .isCharVectorOrNull(value, "color_title")
    if (!is.null(value)) {
      .stopifnot(
        length(value) == 1 ||
          length(value) == length(color_nodes_by),
        "color_title", "must have length 1 or the same length as color_nodes_by"
      )
    }
  }
}

.checkSizeNodesBy <- function(value, data) {
  .isCharOrNumericOrNull(value, "size_nodes_by")
  if (!is.null(value)) {
    .hasLength1(value, "size_nodes_by")
    if (is.character(value)) {
      .isDataCol(data, value, "size_nodes_by")
    } else if (is.numeric(value)) {
      .isPos(value, "size_nodes_by")
      if (length(value) == 2) {
        .stopifnot(
          value[[1]] <= value[[2]],
          "first entry of size_nodes_by", "cannot be greater than the second entry"
        )
      }
    }
  }
}

.checkNodeSizeLimits <- function(value) {
  .isNumericOrNull(value, "node_size_limits")
  if (!is.null(value)) {
    .hasLength2(value, "node_size_limits")
    .isPos(value[[1]], "values for node_size_limits")
    .isPos(value[[2]], "values for node_size_limits")
    .stopifnot(
      value[[1]] <= value[[2]],
      "first entry of node_size_limits", "cannot be greater than the second entry"
    )
  }
}

.checkStatsToInclude <- function(value) {
  .nonNull(value, "stats_to_include")
  if (typeof(value) %in% c("logical", "list")) {
    stopifnot(
      "value for `stats_to_include` does not match required format (see `?chooseNodeStats`)" =
        names(value) == names(chooseNodeStats())
    )
  } else {
    .isString0(value, "stats_to_include")
    stopifnot(
      "`stats_to_include` must be \"all\", \"cluster_id_only\", or generated by `chooseNodeStats`" =
        value %in% c("all", "cluster_id_only")
    )
  }
}


.checkClusterFun <- function(value, argname = "cluster_fun") {
  .stopifnot(
    isTRUE(all.equal(value, cluster_edge_betweenness)) ||
      isTRUE(all.equal(value, cluster_fast_greedy)) ||
      isTRUE(all.equal(value, cluster_fluid_communities)) ||
      isTRUE(all.equal(value, cluster_infomap)) ||
      isTRUE(all.equal(value, cluster_label_prop)) ||
      isTRUE(all.equal(value, cluster_leading_eigen)) ||
      isTRUE(all.equal(value, cluster_leiden)) ||
      isTRUE(all.equal(value, cluster_louvain)) ||
      isTRUE(all.equal(value, cluster_optimal)) ||
      isTRUE(all.equal(value, cluster_spinglass)) ||
      isTRUE(all.equal(value, cluster_walktrap)),
    "cluster_fun",
    "must be a valid clustering algorithm. See ?clustering_algorithms"

  )

}


.checkIgraphAgainstData <- function(igraph, data) {
  .stopifnot(
    length(igraph) == nrow(data),
    "", "number of nodes in `igraph` does not match number of rows in `data`"
  )
}

.checkIgraphAgainstMatrix <- function(igraph, mat) {
  .stopifnot(
    length(igraph) == nrow(mat),
    "", "number of nodes in `igraph` does not match dimensions of `adjacency_matrix`"
  )
}
.checkDataAgainstMatrix <- function(data, mat) {
  .stopifnot(
    nrow(data) == nrow(mat),
    "", "number of data rows does not match dimensions of the adjacency matrix"
  )
}