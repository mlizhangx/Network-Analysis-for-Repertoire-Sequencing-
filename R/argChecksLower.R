
# Universal Checks ------------------------------------------------------------

.stopifnot <- function(condition, argname, message) {
  if (!condition) { stop(paste(argname, message)) }
}

.nonNull <- function(value, argname) {
  if (is.null(value)) { stop(paste(argname, "is required but value is NULL")) }
}

.noNAs <- function(value, argname) {
  if (sum(is.na(value)) > 0) { stop(paste(argname, "contains NAs")) }
}

.hasLength1 <- function(value, argname) {
  if (length(value) != 1) { stop(paste(argname, "must have length 1")) }
}

.hasPosLength <- function(value, argname) {
  if (length(value) == 0) { stop(paste(argname, "must have positive length")) }
}



# Type Checks -------------------------------------------------------------

# checks that an argument is strictly TRUE or FALSE
.isTF <- function(value, argname) {
  .hasLength1(value, argname)
  .stopifnot(sum(value %in% c(TRUE, FALSE)) == 1,
             argname, "must evaluate to \"TRUE\" or \"FALSE\"")
}


# is a string
.isString0 <- function(value, argname) {
  .hasLength1(value, argname)
  .stopifnot(is.character(value), argname, "must be of type 'character'")
}
# is a string and is not NULL
.isString <- function(value, argname) {
  .nonNull(value, argname)
  .isString0(value, argname)
}
# is a string or is NULL
.isStringOrNull <- function(value, argname) {
  if (!is.null(value)) { .isString0(value, argname) }
}
# check string or expression
.isStringOrExpr0 <- function(value, argname) {
  .hasLength1(value, argname)
  .stopifnot(is.character(value) | is.expression(value),
             argname, "must be of type 'character' or 'expression")
}
.isStringOrExpr <- function(value, argname) {
  .nonNull(value, argname)
  .isStringOrExpr0(value, argname)
}
# check string, expression or NULL
.isStringExprOrNull <- .isStringOrExprOrNull <- function(value, argname) {
  if (!is.null(value)) { .isStringOrExpr0(value, argname) }
}


# checks that argument is nonnegative and finite
.isNonneg0 <- function(value, argname) {
  .hasLength1(value, argname)
  .stopifnot(is.numeric(value), argname, "must be of type 'numeric'")
  .stopifnot(is.finite(value), argname, "must be finite")
  .stopifnot(value >= 0, argname, "must be nonnegative")
}
.isNonneg <- function(value, argname) {
  .nonNull(value, argname)
  .isNonneg0(value, argname)
}
.isNonnegOrNull <- function(value, argname) {
  if (!is.null(value)) { .isNonneg0(value, argname) }
}

# checks that argument is strictly positive and finite
.isPos0 <- function(value, argname) {
  .hasLength1(value, argname)
  .stopifnot(is.numeric(value), argname, "must be of type 'numeric'")
  .stopifnot(is.finite(value), argname, "must be finite")
  .stopifnot(value > 0, argname, "must be strictly positive")
}
.isPos <- function(value, argname) {
  .nonNull(value, argname)
  .isPos0(value, argname)
}
.isPosOrNull <- function(value, argname) {
  if (!is.null(value)) { .isPos0(value, argname) }
}

# checks that argument is a strictly positive integer
.isPosInt0 <- function(value, argname) {
  .isPos(value, argname)
  .stopifnot(value %% 1 == 0, argname, "must be an integer")
}
.isPosInt <- function(value, argname) {
  .nonNull(value, argname)
  .isPosInt0(value, argname)
}
.isPosIntOrNull <- function(value, argname) {
  if (!is.null(value)) { .isPosInt0(value, argname) }
}

# check of type character
.isChar <- function(value, argname) {
  .stopifnot(is.character(value), argname, "is of a type other than 'character'")
}

# check of type character or numeric
.isCharOrNumeric <- function(value, argname) {
  .stopifnot(is.character(value) | is.numeric(value), argname, "is of a type other than 'character' or 'numeric'")
}


# Custom Checks -----------------------------------------------------------



.isRepSeqData <- function(data) {
  # has at least two rows
  stopifnot("need at least two data rows" = nrow(data) > 1)
}

.isDistType <- function(value) {
  .isString(value, "dist_type")
  stopifnot("Invalid option for `dist_type` argument" =
              value %in% c("levenshtein", "Levenshtein, lev, Lev, l, L",
                           "hamming", "Hamming", "ham", "Ham", "h", "H"))
}


# type = "network" for SaveNetwork()
# type = "generic" for .saveDataGeneric(), used in findAssociatedClones()
.isOutputType <- function(value, type = "network") {
  .isString(value, "output_type")
  valid_types <- switch(type,
                        "network" = c("individual", "rds", "rda"),
                        "generic" = c("rds", "rda", "csv", "tsv", "table"))
  if (!value %in% valid_types) {
    warning("invalid option for `output_type`, defaulting to 'rda'")
  }
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
    stop(paste(argname, "value has type other than 'character' or 'numeric'"))
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
      for (i in length(colref)) { .isDataCol0(data, colref[[i]], argname) }
    }
  }
}

# Checks that column specified for sequence column is valid
.isSeqCol <- function(data, colref) {
  .nonNull(colref, "`seq_col`")
  if (length(colref) == 1) {
    .isDataCol0(data, colref, "`seq_col`")
  } else if (length(colref) == 2) {
    .isDataCol0(data, colref[[1]], "first entry of `seq_col`")
    .isDataCol0(data, colref[[2]], "second entry of `seq_col`")
  } else {
    stop("`seq_col` must have length 1 or 2")
  }

}



# checks that value is NULL, "auto", or valid column reference
.checkColorNodesBySingle <- function(value, data, node_stats, cluster_stats) {
  .hasLength1(value, "each value in `color_nodes_by`")
  if (is.character(value)) {
    valid_names <- names(data)
    if (node_stats) {
      valid_names <- c(valid_names, names(chooseNodeStats()))
      valid_names <- valid_names[-length(valid_names)]  # remove "all_stats"
    }
    if (cluster_stats) {
      valid_names <-
        c(valid_names,
          "cluster_id", "node_count",
          "eigen_centrality_eigenvalue", "eigen_centrality_index",
          "closeness_centrality_index", "degree_centrality_index",
          "edge_density", "assortativity", "global_transitivity",
          "diameter_length", "seq_w_max_count", "max_count", "agg_count",
          "seq_w_max_degree", "max_degree", "mean_degree", "mean_seq_length")
    }
    stopifnot(
      "`color_nodes_by` specifies one or more column names not present in `data`" =
        value %in% valid_names)
  } else if (is.numeric(value)) {
    stopifnot(
      "`color_nodes_by` specifies one or more column numbers outside of the valid range" =
        value >= 1 && value <= ncol(data))
  }
}

.checkColorNodesBy <- function(value, data, node_stats, cluster_stats) {
  if (!is.null(value)) {
    .hasPosLength(value, "color_nodes_by")
    .isCharOrNumeric(value, "color_nodes_by")
    if (length(value) == 1) {
      if (value != "auto") {
        .checkColorNodesBySingle(value, data, node_stats, cluster_stats)
      }
    } else { # multiple values; check each individually
      for (i in seq_along(value)) {
        .checkColorNodesBySingle(value[[i]], data, node_stats, cluster_stats)
      }
    }
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