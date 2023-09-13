# NAIR: Network Analysis of Immune Repertoire
# Copyright (C) 2023 Li Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Conventions -------------------------------------------------------------

# Checks beginning with `.is` return TRUE or FALSE
# Checks beginning with `.MUST` raise an error or return NULL (hard checks)
# Checks beginning with `.check` return the main argument or a default value
#                                                               (soft checks)



# String Content ----------------------------------------------------------

.isValidFilenamePart <- function(x) {
  # contains only letters, numbers, hyphens and underscores;
  # does not begin or end with a hyphen or underscore
  isTRUE(nchar(x) > 0) &&
  isFALSE(grepl("\\W", gsub("-", "", x))) &&
    isFALSE(grepl("^[-_]", x)) && isFALSE(grepl("[-_]$", x))
}

.isValidFilename <- function(x) {
  # contains only letters, numbers, hyphens, underscores and at most one period;
  # does not begin or end with a hyphen, underscore or period
  isTRUE(nchar(x) > 0) &&
  isFALSE(grepl("\\W", gsub("[.-]", "", x))) &&
    isFALSE(lengths(regmatches(x, gregexpr("\\.", x))) > 1) &&
    isFALSE(grepl("^[.-_]", x)) && isFALSE(grepl("[.-_]$", x))
}

.sanitizeFilenamePart <- function(x) {
  # Replace illegal characters with underscores;
  # remove any leading or trailing hyphens or underscores
  pieces <- unlist(strsplit(x, split = "-"))
  pieces <- gsub("\\W", "_", pieces)
  pieces <- pieces[!pieces %in% ""]
  pieces <- pieces[!grepl("^_+$", pieces)]
  x <- paste(pieces, collapse = "-")
  while (grepl("^_", x)) { x <- substr(x, 2, nchar(x)) }
  while (grepl("_$", x)) { x <- substr(x, 1, nchar(x) - 1) }
  x
}


# Lists and Data Frames ----------------------------------------------------


.isNamedList <- function(x) { is.list(x) && !is.null(names(x)) }
.MUST.isNamedList <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isNamedList(x), name, "must be a named list")
}

.MUST.isDataFrame <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(is.data.frame(x), name, "could not be coerced to a data frame")
}

.MUST.hasMultipleRows <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(nrow(x) > 1, name, "must contain at least two rows")
}


# Column References -------------------------------------------------------

.isDataColref <- function(x, data) {
  (.isString(x) && x %in% names(data)) ||
    (.isPosInt(x) && x <= ncol(data))

}
.MUST.isDataColref <- function(x, data, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isDataColref(x, data), name, "must specify a data column")
}

.isDataColrefs <- function(x, data) {
  (.hasLength1(x) && .isDataColref(x, data)) ||
    (length(x) > 1 && is.vector(x) &&
       all(sapply(x, .isDataColref, data = data))
     )
}
.MUST.isDataColrefs <- function(x, data, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isDataColrefs(x, data), name,
             "must specify one or more data columns"
  )
}
.checkDataColrefs <- function(x, data, default = NULL, ornull = TRUE) {
  # Replace column references with default value if not all valid
  if ((!ornull && is.null(x)) || !is.null(x) && !.isDataColrefs(x, data)) {
    warning(
      "one or more values in ", sQuote(deparse(substitute(x))),
      " are not valid column references for ",
      sQuote(deparse(substitute(data))), ". ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}

.checkDataForColref <- function(x, data, check, check2 = NULL, default = NULL) {
  # Ensure specified data column passes given check(s)
  is_col <- .isDataColref(x, data)
  invalid <- !is_col || (is_col && !check(data[[x]]))
  if (!is.null(check2)) {
    invalid <- invalid || (is_col && !check2(data[[x]]))
  }
  if (invalid) {
    warning(
      "value for ", sQuote(deparse(substitute(x))),
      " does not reference a column of ",
      sQuote(deparse(substitute(data))), " with the correct properties. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}




# AIRR-Seq Data -----------------------------------------------------------


.isValidSeqVector <- function(x) {
  x <- as.vector(x, mode = "character")
  is.vector(x, mode = "character") && .hasPosLength(x) && sum(!is.na(x)) > 0
}
.MUST.isValidSeqVector <- function(x, name = NULL, sep = " ") {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isValidSeqVector(x),
    name,
    "must be coercible to a character vector with at least one non-NA value",
    sep = sep
  )
}

.isSeqColref <- function(x, data) {
  .isDataColref(x, data) && .isValidSeqVector(data[[x]])
}
.MUST.isSeqColref <- function(x, data, name = NULL, data_name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(data_name)) { name <- deparse(substitute(data)) }
  .stopifnot(.isSeqColref(x, data),
             name, "does not reference a column of", sQuote(data_name),
             "that is coercible to a",
             "character vector with at least one non-NA value"
  )
}

.isSeqColrefs <- function(x, data) {
  length(x) %in% c(1, 2) && .isDataColrefs(x, data) &&
    all(sapply(x, function(x) { .isValidSeqVector(data[[x]]) }))
}
.MUST.isSeqColrefs <- function(x, data, name = NULL, data_name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(data_name)) { name <- deparse(substitute(data)) }
  .stopifnot(
    .isSeqColrefs(x, data),
    name, "does not reference one or two columns of", sQuote(data_name),
    "that are each coercible to a",
    "character vector with at least one non-NA value"
  )
}

.isCountColref <- function(x, data) {
  .isDataColref(x, data) && is.numeric(data[[x]])
}
.MUST.isCountColref <- function(x, data, name = NULL, data_name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(data_name)) { name <- deparse(substitute(data)) }
  .stopifnot(.isCountColref(x, data),
             name, "does not reference a column of", sQuote(data_name),
             "that contains numeric values"
  )
}
.checkCountCol <- function(x, data, default = NULL) {
  if (!is.null(x) && !.isCountColref(x, data)) {
    warning(
      "value for ", sQuote("count_col"), " is not a valid column reference ",
      "or does not specify a numeric variable of ",
      sQuote(deparse(substitute(data))), ". ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}

# Network Objects ----------------------------------------------------------



.isGgraph <- function(x) { inherits(x, "ggraph") }
.MUST.isGgraph <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isGgraph(x), name, "must be of class", dQuote("ggraph"))
}

.isIgraph <- function(x) { inherits(x, "igraph") }
.MUST.isIgraph <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isIgraph(x), name, "must be of class", dQuote("igraph"))
}

.isAdjacencyMatrix <- function(x) {
  .hasPosLength(x) &&
    (is.matrix(x) || inherits(x, "sparseMatrix")) &&
    isTRUE(dim(x)[[1]] == dim(x)[[2]]) &&
    Matrix::isSymmetric(x, check.attributes = FALSE) &&
    all(unique(as.vector(x)) %in% c(0, 1))
}
.MUST.isAdjacencyMatrix <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isAdjacencyMatrix(x),
             name, "must be a symmetric matrix",
             "(including sparse formats from the Matrix package)",
             "containing only ones and zeros"
  )
}

.isLayout <- function(x) {
  is.matrix(x) && is.numeric(x) && ncol(x) == 2 && nrow(x) > 0
}
.MUST.isLayout <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isLayout(x), name, "must be a nonempty two-column numeric matrix")
}

.isPlotlist <- function(x) {
  out <- .isNamedList(x)
  if (.hasElem(x, "graph_layout")) {
    non_layout_indices <- which(names(x) != "graph_layout")
    out <- out && .hasPosLength(non_layout_indices) && .isLayout(x$graph_layout)
  } else {
    non_layout_indices <- 1:length(x)
  }
  out && all(sapply(x[non_layout_indices], .isGgraph))
}
.MUST.isPlotlist <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isPlotlist(x),
    name, "must be a named list containing objects of class", dQuote("ggraph"),
    "and possibly a two-column numeric matrix named", dQuote("graph_layout")
  )
}


# Correspondence of Network Objects -----------------------------------------------



.doesIgraphMatchData <- function(x, data) { length(x) == nrow(data) }
.MUST.doesIgraphMatchData <- function(x, data, name = NULL, data_name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(data_name)) { name <- deparse(substitute(data)) }
  .stopifnot(
    .doesIgraphMatchData(x, data),
    name, "must have length equal to row dimension of", sQuote(data_name)
  )
}

.doesIgraphMatchMatrix <- function(x, mat) { length(x) == nrow(mat) }
.MUST.doesIgraphMatchMatrix <- function(x, mat, name = NULL, matname = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(matname)) { name <- deparse(substitute(mat)) }
  .stopifnot(
    .doesIgraphMatchMatrix(x, mat),
    name, "must have length equal to row dimension of", sQuote(matname)
  )
}

.doesDataMatchMatrix <- function(x, mat) { nrow(x) == nrow(mat) }
.MUST.doesDataMatchMatrix <- function(x, mat, name = NULL, matname = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(matname)) { name <- deparse(substitute(mat)) }
  .stopifnot(.doesDataMatchMatrix(x, mat),
             name, "and", sQuote(matname), "must have the same row dimension"
  )
}

.doesPlotMatchData <- function(x, data) {
  # x is a ggraph, data is a data frame
  nrow(extractLayout(x)) == nrow(data)
}
.MUST.doesPlotMatchData <- function(x, data, name = NULL, data_name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(data_name)) { name <- deparse(substitute(data)) }
  .stopifnot(
    .doesPlotMatchData(x, data),
    name, "has node count not equal to the row dimension of", sQuote(data_name)
  )
}
.doesLayoutMatchData <- function(x, data) {
  .isLayout(x) && nrow(x) == nrow(data)
}
.MUST.doesLayoutMatchData <- function(x, data, name = NULL, data_name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  if (is.null(data_name)) { name <- deparse(substitute(data)) }
  .stopifnot(.doesLayoutMatchData(x, data),
             name, "must be NULL or a two-column matrix",
             "with the same row dimension as", sQuote(data_name)
  )
}

# Network List ------------------------------------------------------------


.hasIgraph <- function(x) {
  .isNamedList(x) && .hasElem(x, "igraph") && .isIgraph(x$igraph)
}
.MUST.hasIgraph <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.hasIgraph(x),
             name, "must contain an object of class", dQuote("igraph"),
             "named", dQuote("igraph")
  )
}

.hasAdjacencyMatrix <- function(x) {
  .isNamedList(x) && .hasElem(x, "adjacency_matrix") &&
    .isAdjacencyMatrix(x$adjacency_matrix)
}
.MUST.hasAdjacencyMatrix <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .hasAdjacencyMatrix(x),
    name,
    "must contain a symmetric matrix named", dQuote("adjacency_matrix"),
    "(including sparse formats from the Matrix package)",
    "that contains only ones and zeros"
  )
}


.hasNodeData <- function(x) {
  .isNamedList(x) && .hasElem(x, "node_data") &&
    .hasPosLength(x$node_data) && is.data.frame(x$node_data)
}
.MUST.hasNodeData <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .hasNodeData(x),
    name, "must contain a nonempty data frame named", dQuote("node_data")
  )
}

.hasDetails <- function(x) {
  .isNamedList(x) && .hasElem(x, "details") &&
    .hasPosLength(x$details) && .isNamedList(x$details)
}
.MUST.hasDetails <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .hasDetails(x),
    name, "must contain a nonempty named list named", dQuote("details")
  )
}

.hasClusterData <- function(x) {
  .isNamedList(x) && .hasElem(x, "cluster_data") &&
    .hasPosLength(x$cluster_data) && is.data.frame(x$cluster_data) &&
    .hasElem(x$cluster_data, "cluster_id") &&
    .isIntegerVector(x$cluster_data$cluster_id, factor_ok = TRUE)
}
.MUST.hasClusterData <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .hasClusterData(x),
    name, "must contain a nonempty data frame named", dQuote("cluster_data"),
    "with an integer-valued variable named", dQuote("cluster_id")
  )
}


.hasPlots <- function(x) {
  .isNamedList(x) && .hasElem(x, "plots") && .isPlotlist(x$plots)
}
.MUST.hasPlots <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .hasPlots(x),
    name, "must contain a named list named", dQuote("plots"),
    "containing objects of class", dQuote("ggraph"),
    "and possibly a two-column numeric matrix named", dQuote("graph_layout")
  )
}


.isBaseNetworkOutput <- function(x) {
  .hasIgraph(x) && .hasAdjacencyMatrix(x) && .hasNodeData(x) &&
    .doesIgraphMatchData(x$igraph, x$node_data) &&
    .doesIgraphMatchMatrix(x$igraph, x$adjacency_matrix) &&
    .doesDataMatchMatrix(x$node_data, x$adjacency_matrix)
}
.MUST.isBaseNetworkOutput <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isBaseNetworkOutput(x),
    name, "must be a named list containing elements",
    paste0(dQuote("igraph"), ","),
    dQuote("adjacency_matrix"), "and", paste0(dQuote("node_data"), ","),
    "all corresponding to the same network"
  )
}




