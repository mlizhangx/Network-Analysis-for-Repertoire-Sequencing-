
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
    stop(paste(argname, "must not contain NA/NaNs"))
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

# combine with other checks to allow NULL argument values
# only works with checks that have only value, argname as inputs
.orNull <- function(check, value, argname) {
  if (!is.null(value)) {
    check(value, argname)
  }
}

.isLogical <- function(value, argname) {
  .stopifnot(
    is.logical(value),
    argname, "must be of type logical"
  )
  .noNAs(value, argname)
}

.isChar <- function(value, argname) {
  .stopifnot(
    is.character(value), argname,
    "must be of type character"
  )
}

.isNumeric <- function(value, argname) {
  .stopifnot(
    is.numeric(value),
    argname, "must be of type numeric"
  )
}

.isCharOrNumeric <- function(value, argname) {
  .stopifnot(
    is.character(value) || is.numeric(value),
    argname, "must be of type character or numeric"
  )
}

.isCharOrLogical <- function(value, argname) {
  .stopifnot(
    is.character(value) || is.logical(value),
    argname, "must be of type character or logical"
  )
  .noNAs(value, argname)
}

.isTF <- function(value, argname) {
  .isLogical(value, argname)
  .hasLength1(value, argname)
  .stopifnot(sum(value %in% c(TRUE, FALSE)) == 1,
             argname, "must evaluate to \"TRUE\" or \"FALSE\""
  )
}
.isTFOrAuto <- function(value, argname) {
  .nonNull(value, argname)
  .isCharOrLogical(value, argname)
  .hasLength1(value, argname)
  if (is.character(value)) {
    .stopifnot(
      value == "auto",
      argname, "must be TRUE, FALSE, or \"auto\""
    )
  } else {
    .stopifnot(
      sum(value %in% c(TRUE, FALSE)) == 1,
      argname, "must be TRUE, FALSE, or \"auto\""
    )
  }
}

.isFinite <- function(value, argname) {
  .isNumeric(value, argname)
  .hasPosLength(value, argname)
  .stopifnot(
    all(is.finite(value)),
    argname, "must contain finite values"
  )
}

.isNonneg <- function(value, argname) {
  .nonNull(value, argname)
  .isNumeric(value, argname)
  .hasLength1(value, argname)
  .stopifnot(
    all(is.finite(value)),
    argname, "must contain finite values"
  )
  .stopifnot(
    value >= 0,
    argname, "must be nonnegative"
  )
}

.isPos <- function(value, argname) {
  .nonNull(value, argname)
  .isNumeric(value, argname)
  .hasLength1(value, argname)
  .isFinite(value, argname)
  .stopifnot(
    value > 0,
    argname, "must be strictly positive"
  )
}

.isInt <- function(value, argname) {
  .nonNull(value, argname)
  .isNumeric(value, argname)
  .hasLength1(value, argname)
  .isFinite(value, argname)
  .stopifnot(value %% 1 == 0,
             argname, "must be integer-valued"
  )
}

.isPosInt <- function(value, argname) {
  .nonNull(value, argname)
  .isNumeric(value, argname)
  .hasLength1(value, argname)
  .isFinite(value, argname)
  .stopifnot(value %% 1 == 0,
             argname, "must be integer-valued"
  )
  .stopifnot(
    value > 0,
    argname, "must be strictly positive"
  )
}

.isString <- function(value, argname) {
  .nonNull(value, argname)
  .isChar(value, argname)
  .hasLength1(value, argname)
}

.isCharOrNumericScalar <- function(value, argname) {
  .nonNull(value, argname)
  .isCharOrNumeric(value, argname)
  .hasLength1(value, argname)
  if (is.numeric(value)) {
    .isFinite(value, argname)
  }
}

.isCharVector <- function(value, argname) {
  .nonNull(value, argname)
  .stopifnot(
    is.vector(value, mode = "character"),
    argname, "must be a character vector"
  )
  .hasPosLength(value, argname)
  .noNAs(value, argname)
}

.isCharOrNumericVector <- function(value, argname) {
  .nonNull(value, argname)
  .stopifnot(
    is.vector(value, mode = "character") || is.vector(value, mode = "numeric"),
    argname, "must be a character or numeric vector"
  )
  .hasPosLength(value, argname)
  .noNAs(value, argname)
  if (is.numeric(value)) {
    .isFinite(value, argname)
  }
}


# Data Frames and Lists ----------------------------------------------------


.hasAtLeastTwoRows <- function(data) {
  stopifnot("need at least two data rows" = nrow(data) > 1)
}

.isDataCol <- function(data, colref, argname) {
  .nonNull(colref, argname)
  .isCharOrNumeric(colref, argname)
  .hasLength1(colref, argname)
  if (is.character(colref)) {
    .stopifnot(
      colref %in% names(data),
      argname, "must specify a valid column within the data")
  } else if (is.numeric(colref)) {
    .isInt(colref, paste(argname, "is of type numeric and hence"))
    .stopifnot(
      colref >= 1 && colref <= ncol(data),
      argname, "must specify a valid column within the data")
  }
}
.isDataColOrNull <- function(data, colref, argname) {
  if (!is.null(colref)) {
    .isDataCol(data, colref, argname)
  }
}

.isDataCols <- function(data, colref, argname) {
  .nonNull(colref, argname)
  .isCharOrNumeric(colref, argname)
  .hasPosLength(colref, argname)
  if (length(colref) == 1) {
    .isDataCol(data, colref, argname)
  } else if (length(colref) > 1) {
    .noNAs(colref, argname)
    for (i in length(colref)) {
      .isDataCol(
        data, colref[[i]],
        paste("each element of", argname)
      )
    }
  }
}
.isDataColsOrNull <- function(data, colref, argname) {
  if (!is.null(colref)) {
    .isDataCols(data, colref, argname)
  }
}

.isValidSeqVector <- function(value) {
  value <- as.vector(value, mode = "character")
  if (!is.vector(value, mode = "character")) {
    stop("specified column or vector of receptor sequences must be coercible to a character vector")
  }
  .hasPosLength(value, "specified column or vector of receptor sequences")
  .stopifnot(
    sum(is.na(value)) < length(value),
    "specified column or vector of receptor sequences",
    "contains only NA values after being coerced to a character vector"
  )
}

.isSeqCol <- function(data, colref) {
  .nonNull(colref, "seq_col")
  .stopifnot(
    length(colref) %in% c(1, 2),
    "seq_col", "must have length 1 or 2"
  )
  .isDataCols(data, colref, "seq_col")
  if (length(colref) == 1) {
    .isValidSeqVector(data[[colref]])
  } else if (length(colref) == 2) {
    .isValidSeqVector(data[[colref[[1]]]])
    .isValidSeqVector(data[[colref[[2]]]])
  }
}

.hasElement <- function(value, argname, element_name) {
  .stopifnot(
    element_name %in% names(value),
    argname, paste("does not contain an element named", element_name)
  )
}

# Network Output ----------------------------------------------------------

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
    argname, "contains values other than 0 or 1"
  )
}

.hasNodeAndClusterData <- function(value, argname) {
  .isList(value, argname)
  .hasPosLength(value, argname)
  .hasElement(value, argname, "node_data")
  .hasElement(value, argname, "cluster_data")
  .isDataFrame(value$node_data, paste0(argname, "$node_data"))
  .isDataFrame(value$cluster_data, paste0(argname, "$cluster_data"))
}

.checkIgraphAgainstData <- function(igraph, data) {
  .stopifnot(
    length(igraph) == nrow(data),
    "", "number of nodes in igraph does not match number of rows in data"
  )
}

.checkIgraphAgainstMatrix <- function(igraph, mat) {
  .stopifnot(
    length(igraph) == nrow(mat),
    "", "number of nodes in igraph does not match dimensions of adjacency matrix"
  )
}
.checkDataAgainstMatrix <- function(data, mat) {
  .stopifnot(
    nrow(data) == nrow(mat),
    "", "number of data rows does not match dimensions of the adjacency matrix"
  )
}

.isPlotlist <- function(value, argname) {
  .isList(value, argname)
  .hasPosLength(value, argname)
  for (i in 1:length(value)) {
    .isGgraph(value[[i]], paste("each element of", argname))
  }
}

.isBaseNetworkOutput <- function(value, argname) {
  .isList(value, argname)
  .hasPosLength(value, argname)
  .hasElement(value, argname, "node_data")
  .hasElement(value, argname, "igraph")
  .hasElement(value, argname, "adjacency_matrix")
  .isDataFrame(value$node_data, paste0(argname, "$node_data"))
  .isIgraph(value$igraph, paste0(argname, "$igraph"))
  .isAdjacencyMatrix(
    value$adjacency_matrix, paste0(argname, "$adjacency_matrix")
  )
  .checkIgraphAgainstData(value$igraph, value$node_data)
  .checkIgraphAgainstMatrix(value$igraph, value$adjacency_matrix)
}

# Specific Arguments -----------------------------------------------------------

.isDistType <- function(value) {
  .isString(value, "dist_type")
  stopifnot(
    "Invalid option for dist_type argument" =
      value %in% c(
        "levenshtein", "Levenshtein", "lev", "Lev", "l", "L",
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
# type = "generic" for .saveDataGeneric()
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

.checkColorNodesBySingle <- function(value, data, node_stats) {
  valid_names <- names(data)
  if (node_stats) {
    valid_names <- c(valid_names, names(chooseNodeStats()))
    valid_names <- valid_names[-length(valid_names)]  # remove "all_stats"
  }
  stopifnot(
    "color_nodes_by specifies one or more variables not present in data or among the node-level network properties to be computed" =
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
  .stopifnot(
    value %in% valid_names,
    "color_scheme", "contains one or more values which are not supported. See help file for plotNetworkGraph()"
  )
}

.checkColorTitle <- function(value, color_nodes_by) {
  if (!is.null(color_nodes_by)) {
    .orNull(.isCharVector, value, "color_title")
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
  .orNull(.isCharOrNumeric, value, "size_nodes_by")
  if (!is.null(value)) {
    .hasLength1(value, "size_nodes_by is non-null and hence")
    if (is.character(value)) {
      .isDataCol(data, value, "size_nodes_by is of type character and hence")
    } else if (is.numeric(value)) {
      .isPos(value, "size_nodes_by is of type numeric and hence")
    }
  }
}

.checkNodeSizeLimits <- function(value) {
  .orNull(.isNumeric, value, "node_size_limits")
  if (!is.null(value)) {
    .hasLength2(value, "node_size_limits is non-null and hence")
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
    .stopifnot(
      length(names(value)) == length(names(chooseNodeStats())) &&
        all(names(value) == names(chooseNodeStats())),
      "value for stats_to_include",
      "does not match required format. See help file for chooseNodeStats()"
    )
  } else {
    .isString(value, "stats_to_include")
    .stopifnot(
      value %in% c("all", "cluster_id_only"),
      "value for stats_to_include",
      "does not match required format. See help file for chooseNodeStats()"
    )
  }
}

.checkClusterFun <- function(value, argname = "cluster_fun") {
  .nonNull(value, argname)
  .hasLength1(value, argname)
  if (is.character(value)) {
    .stopifnot(
      value %in% c(
        "cluster_edge_betweenness", "cluster_fast_greedy",
        "cluster_fluid_communities", "cluster_infomap",
        "cluster_label_prop", "cluster_leading_eigen",
        "cluster_leiden", "cluster_louvain",
        "cluster_optimal", "cluster_spinglass", "cluster_walktrap"
      ),
      argname, "must be a valid clustering algorithm. See help topic \"clustering_algorithms\""
    )
  } else {
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
      argname, "must be a valid clustering algorithm. See help topic \"clustering_algorithms\""
    )
  }
}
