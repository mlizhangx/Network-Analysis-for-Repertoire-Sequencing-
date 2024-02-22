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


# File Input Arguments ----------------------------------------------------


.isInputType <- function(x) {
  choices <- c("csv", "csv2", "table", "tsv", "txt", "rds", "rda")
  .isString(x) && pmatch(x, choices, nomatch = 0)
}
.MUST.isInputType <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isInputType(x),
             name, "must be one of:",
             paste(dQuote("rds"), dQuote("rda"), dQuote("csv"), dQuote("csv2"),
                   dQuote("tsv"), dQuote("table"), sep = ", ")
  )
}

.checkReadArgs <- function(x, header, sep) {
  if (!is.null(x)) {
    choices <- names(formals(utils::read.table))
    choices <- choices[2:length(choices)]
    matches <- rep(TRUE, length(x))
    for (i in 1:length(x)) {
      if (!names(x)[i] %in% choices) {
        matches[i] <- FALSE
        warning("dropping argument ", sQuote(names(x)[i]),
                " from ", sQuote("read.args"), " as it is not an optional ",
                "argument to ", sQuote("utils::read.table()")
        )
      }
      x <- x[matches]
    }
    if (!"header" %in% names(x)) { x$header <- header }
    if (!"sep" %in% names(x)) { x$sep <- sep }
  }
  x
}

.checkargs.InputFiles <- function(
    file_list, input_type, data_symbols, header, sep, read.args = NULL
) {
  .stopifnot(.isCharVector(file_list) || is.list(file_list),
             "file_list",
             "must be a character vector",
             "or a list of character strings and connections"
  )
  if (is.list(file_list)) {
    .stopifnot(
      all(sapply(file_list, inherits, what = c("connection", "character"))),
      "file_list",
      "contains elements other than character strings and connections"
    )
    .stopifnot(
      all(sapply(file_list, .hasLength1)),
      "file_list",
      "contains elements other than character strings and connections"
    )
    string_positions <- which(sapply(file_list, inherits, what = "character"))
  } else {
    string_positions <- 1:length(file_list)
  }
  .stopifnot(isTRUE(all(sapply(file_list[string_positions], file.exists))),
             "file_list", "specifies one or more nonexistent files"
  )
  .stopifnot(length(file_list) == length(unique(file_list)),
             "file_list", "contains duplicate values"
  )
  .MUST.isInputType(input_type)
  if (input_type == "rda") {
    .MUST.isCharVector(data_symbols, "data_symbols")
    .stopifnot(
      length(data_symbols) %in% c(1, length(file_list)),
      "data_symbols", "must have length 1 or equal to that of",
      sQuote("file_list")
    )
  }
  if (input_type %in% c("csv", "csv2", "table", "tsv", "txt")) {
    .MUST.isTF(header)
    .MUST.isString(sep)
    if (!is.null(read.args)) {
      .MUST.isNamedList(read.args)
    }
  }
}

.checkIDs <- function(x, len, default = NULL,
                      ornull = FALSE, allow_dupes = FALSE
) {
  if (ornull && is.null(x)) {
    return(x)
  }
  if (!.isCharOrNumericVector(x, factor_ok = TRUE) || length(x) != len ||
      (!allow_dupes && sum(duplicated(x)) > 0)
  ) {
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Using default value instead"
    )
    return(default)
  }
  x
}



# File Output Arguments ---------------------------------------------------

.checkOutputDir <- function(x, name = "output_dir") {
  if (!is.null(x) && !isTRUE(dir.exists(x))) {
    warning("directory ", dQuote(x), " specified for ", sQuote(name),
            " does not exist and could not be created. Output will not be saved"
    )
    return(NULL)
  }
  x
}

.requireOutputDir <- function(x, name = "output_dir") {
  .stopifnot(
    isTRUE(dir.exists(x)),
    NULL, "directory", dQuote(x), "specified for", sQuote(name),
    "does not exist and could not be created"
  )
}

# type = "network" for SaveNetwork()
# type = "generic" for .saveDataGeneric()
.isOutputType <- function(x, type = "network") {
  choices <- switch(
    type,
    "network" = c("individual", "rds", "rda"),
    "findPublicClusters" = c("rds", "rda", "csv"),
    "findAssociatedClones" = c("csv", "tsv", "rds", "rda", "table"),
    "generic" = c("rds", "rda", "csv", "tsv", "table")
  )
  .isString(x) && pmatch(x, choices, nomatch = 0)
}

.checkOutputType <- function(x, type = "network", default = "rds") {
  if (!.isOutputType(x, type)) {
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Defaulting to ", dQuote(default)
    )
    return(default)
  }
  x
}

.checkOutputName <- function(x, default = "MyRepSeqNetwork",
                             argname = "output_name"
) {
  if (!.isNonemptyString(x)) {
    warning(
      "value for ", sQuote(argname), " is not a nonempty character string. ",
      "Defaulting to ", dQuote(default)
    )
    return(default)
  }
  if (!.isValidFilenamePart(x)) {
    x <- .sanitizeFilenamePart(x)
    if (!.isValidFilenamePart(x)) { x <- default }
    warning(
      "value for ", sQuote(argname), " may be unsafe ",
      "for use as a file name prefix. Value changed to ", dQuote(x)
    )
  }
  x
}

.checkOutfileLayout <- function(x, plotlist, default = NULL) {
  if (!is.null(x) && (!.isString(x) || !.hasElem(plotlist, "graph_layout"))) {
    warning(
      sQuote("outfile_layout"), " is non-null but ",
      sQuote("plotlist"), " does not contain a valid layout matrix named ",
      dQuote("graph_layout"), ". No layout will be saved"
    )
    return(default)
  }
  x
}


# Network Analysis Arguments ----------------------------------------------

.isDistType <- function(x) {
  choices <- c("hamming", "levenshtein", "Hamming", "Levenshtein")
  .isString(x) && pmatch(x, choices, nomatch = 0)
}
.MUST.isDistType <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isDistType(x),
             name, "must be", dQuote("hamming"), "or", dQuote("levenshtein"),
             "(abbreviations allowed)"
  )
}

.matchDistType <- function(dist_type) {
  if (pmatch(dist_type, c("levenshtein", "Levenshtein"), 0)) {
    return("levenshtein")
  } else if (pmatch(dist_type, c("hamming", "Hamming"), 0)) {
    return("hamming")
  }
  "hamming"
}

.matchMethod <- function(x, cutoff) {
  if (pmatch(x, "pattern", 0)) {
    if (!any(cutoff == c(0, 1, 2))) {
      warning(
        "pattern algorithm only supports 'dist_cutoff' values of 0, 1 and 2.",
        "Defaulting to ", dQuote("default")
      )
      return("default")
    }
    return("pattern")
  } else if (pmatch(x, "sort", 0)) {
    if (!any(cutoff == c(0, 1))) {
      warning(
        "sort algorithm only supports 'dist_cutoff' value 1.",
        "Defaulting to ", dQuote("default")
      )
      return("default")
    }
    return("sort")
  }
  "default"
}

.checkDistType <- function(x, default = "hamming") {
  if (!.isDistType(x)) {
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  .matchDistType(x)
}

.checkMethod <- function(x, cutoff, default = "default") {
  if (!.isString(x) || !pmatch(x, c("default", "pattern", "sort"), nomatch = 0)) {
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  .matchMethod(x, cutoff)
}

.isStatsToInclude <- function(x) {
  .isLogicalVector(x) && !.hasNAs(x) &&
    isTRUE(all.equal(names(x), names(chooseNodeStats())))
}

.checkStatsToInclude <- function(x, default = chooseNodeStats()) {
  if (.isString(x))  {
    if (x == "all") {
      x <- chooseNodeStats(all_stats = TRUE)
    } else if (x == "cluster_id_only") {
      x <- exclusiveNodeStats(cluster_id = TRUE)
    }
  }
  if (!.isStatsToInclude(x)) {
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Defaulting to ", dQuote(deparse(substitute(default)))
    )
    return(default)
  }
  x
}

.isClusterFun <- function(x) {
  choices <- c(
    "cluster_edge_betweenness", "cluster_fast_greedy",
    "cluster_infomap",
    "cluster_label_prop", "cluster_leading_eigen",
    "cluster_leiden", "cluster_louvain",
    "cluster_optimal", "cluster_spinglass", "cluster_walktrap",
    "edge_betweenness", "betweenness", "fast_greedy", "greedy",
    "infomap",
    "label_prop", "leading_eigen", "prop", "eigen",
    "leiden", "louvain",
    "optimal", "spinglass", "walktrap"
  )
  .isString(x) && pmatch(x, choices, nomatch = 0)
}



# Plotting Arguments ------------------------------------------------------



.checkSizeNodesBy <- function(x, data, default = 0.5) {
  if (!.orNull(.isCharOrNumericScalar, x) ||
      (is.character(x) && !.isDataColref(x, data)) ||
      (is.numeric(x) && !.isPos(x))
  ) {
    warning(
      "value for ", sQuote("size_nodes_by"), " is invalid. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}

.isColorNodesBy <- function(x, data, node_stats = FALSE, cluster_stats = FALSE,
                            plots = TRUE, cluster_id_name = "cluster_id",
                            stats_to_include = chooseNodeStats(),
                            auto_ok = FALSE
) {
  is.null(x) || isFALSE(plots) || (isTRUE(auto_ok) && isTRUE(x == "auto")) ||
    (.isString(x) &&
       .isColorNodesBy1(x, data, node_stats, cluster_stats,
                        cluster_id_name, stats_to_include
       )
    ) ||
    (.isCharVector(x) &&
       all(sapply(x, .isColorNodesBy1, data = data,
                  node_stats = node_stats, cluster_stats = cluster_stats,
                  cluster_id_name = cluster_id_name,
                  stats_to_include = stats_to_include
       ))
    )
}

.isColorNodesBy1 <- function(x, data, node_stats, cluster_stats,
                             cluster_id_name, stats_to_include
) {
  choices <- names(data)
  if (node_stats) {
    if (.isString(stats_to_include)) {
      if (stats_to_include == "all") {
        stats_to_include <- chooseNodeStats(all_stats = TRUE)
      } else if (stats_to_include == "cluster_id_only") {
        stats_to_include <- exclusiveNodeStats(cluster_id = TRUE)
      }
    }
    choices <- c(choices, names(stats_to_include)[stats_to_include])
    if (stats_to_include[["cluster_id"]] && cluster_id_name != "cluster_id") {
      choices <- choices[choices != "cluster_id"]
      choices <- c(choices, cluster_id_name)
    }
  }
  if (cluster_stats) { choices <- c(choices, cluster_id_name) }
  isTRUE(x %in% choices)
}

.checkColorNodesBy <- function(
    x, data, node_stats = FALSE, cluster_stats = FALSE,
    plots = TRUE, cluster_id_name = "cluster_id",
    stats_to_include = chooseNodeStats(), default = NULL, auto_ok = FALSE
) {
  if (!.isColorNodesBy(x, data, node_stats, cluster_stats, plots,
                       cluster_id_name, stats_to_include, auto_ok)
  ) {
    warning(
      "one or more values in ", sQuote("color_nodes_by"), " are invalid. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}

.isColorScheme <- function(x, len) {
  isTRUE(len == 0) ||
    (.isCharVector(x) && length(x) %in% c(1, len) &&
       all(sapply(x, .isColorScheme1))
    )
}
.isColorScheme1 <- function(x) {
  choices <- c(
    "default", "viridis", "magma", "inferno", "plasma", "cividis", "rocket",
    "mako", "turbo", "A", "B", "C", "D", "E", "F", "G", "H"
  )
  choices <- c(choices,
               paste0(choices, "-1"),
               grDevices::hcl.pals()
  )
  isTRUE(x %in% choices)
}

.checkColorScheme <- function(x, color_nodes_by, default = "default") {
  if (!.isColorScheme(x, length(color_nodes_by))) {
    warning("value for ", sQuote("color_scheme"), " is invalid. ",
            "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}


.checkColorTitle <- function(x, color_nodes_by, default = "auto") {
  if (!is.null(color_nodes_by) &&
      (!.orNull(.isCharVector, x) ||
       (!is.null(x) && !length(x) %in% c(1, length(color_nodes_by)))
      )
  ) {
    warning("value for ", sQuote("color_title"), " is invalid. ",
            "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}


.checkNodeSizeLimits <- function(x, default = NULL) {
  if (is.null(x)) {
    return(x)
  }
  if (!(.isNumericVector(x) && .hasLength2(x) && x[[1]] > 0 && x[[2]] >= x[[1]])
  ) {
    warning(
      "value for ", sQuote("node_size_limits"), " is invalid. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}


.checkPlotsAgainstLayout <- function(plots) {
  # plots is a ggraph object
  # assumes plots$graph_layout is a nonempty two-column matrix if it exists
  if (!.hasElem(plots, "graph_layout")) {
    return(NULL)
  }
  non_layout_indices <- which(names(plots) != "graph_layout")
  compare_layout <- function(x) {
    isTRUE(all.equal(plots$graph_layout, as.matrix(x$data[c("x", "y")]),
                     check.attributes = FALSE
    ))
  }
  if (!all(sapply(plots[non_layout_indices], compare_layout))) {
    warning(
      "Graph plot layout contained in ", sQuote("net$plots$graph_layout"),
      " does not match one or more of the plots contained in ",
      sQuote("net$plots"), ". The graph plot layout will be used nonetheless."
    )
  }
  NULL
}


# Associated Clusters -----------------------------------------------------


.checkargs.findAssociatedSeqs <- function(file_list, subject_ids, group_ids,
                                          seq_col
) {
  group_ids <- as.vector(group_ids, mode = "character")
  .orNull(.MUST.isCharOrNumericVector, subject_ids, "subject_ids",
          factor_ok = TRUE
  )
  .MUST.isCharOrNumericVector(group_ids, "group_ids", factor_ok = TRUE)
  .MUST.isCharOrNumericScalar(seq_col, "seq_col")
  .stopifnot(length(file_list) == length(unique(file_list)),
             "file_list", "contains duplicate values"
  )
  .stopifnot(length(file_list) == length(group_ids),
             "file_list", "and", sQuote("group_ids"),
             "have different lengths"
  )
  .stopifnot(length(unique(group_ids)) == 2,
             "group_ids", "must contain exactly two unique values"
  )
  if (!is.null(subject_ids)) {
    .stopifnot(length(file_list) == length(subject_ids),
               "file_list", "and", sQuote("subject_ids"),
               "have different lengths"
    )
  }
}

.checkargs.findAssociatedSeqs2 <- function(
    data, seq_col, sample_col, group_col, data_name, seq_col_name
) {

  .MUST.isDataFrame(data, data_name)
  .MUST.hasMultipleRows(data, data_name)
  .MUST.isSeqColref(seq_col, data, seq_col_name, data_name)
  .MUST.isDataColref(data, group_col, "group_col")
  .MUST.isDataColref(data, sample_col, "sample_col")
  .stopifnot(length(unique(data[[group_col]])) == 2,
             paste0("data[[", group_col, "]]"),
             "must contain exactly two unique values"
  )

}

.checkargs.findAssociatedClones <- function(
    file_list, input_type, data_symbols, header, sep, read.args,
    subject_ids, group_ids, seq_col, assoc_seqs, output_dir
) {

  .checkargs.InputFiles(
    file_list, input_type, data_symbols, header, sep, read.args
  )
  .orNull(.MUST.isCharOrNumericVector, subject_ids, "subject_ids",
          factor_ok = TRUE
  )
  .MUST.isCharOrNumericVector(group_ids, "group_ids", factor_ok = TRUE)
  .MUST.isCharOrNumericScalar(seq_col, "seq_col")
  .MUST.isCharVector(assoc_seqs, "assoc_seqs")
  .stopifnot(length(file_list) == length(unique(file_list)),
             "file_list", "contains duplicate values"
  )
  .stopifnot(length(file_list) == length(group_ids),
             "file_list", "and", sQuote("group_ids"),
             "have different lengths"
  )
  .stopifnot(length(unique(group_ids)) == 2,
             "group_ids", "must contain exactly two unique values"
  )
  if (!is.null(subject_ids)) {
    .stopifnot(length(file_list) == length(subject_ids),
               "file_list", "and", sQuote("subject_ids"),
               "have different lengths"
    )
  }
  .MUST.isString(output_dir, "output_dir")

}


# Deprecated Arguments ----------------------------------------------------
.checkDeprecated.buildPublicClusterNetwork <- function(
    node_stats, stats_to_include, cluster_stats,
    env = caller_env(),
    user_env = caller_env(2)
) {
  if (lifecycle::is_present(node_stats)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "buildPublicClusterNetwork(node_stats)",
      details =
        "All node-level network properties are now automatically computed.",
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(stats_to_include)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "buildPublicClusterNetwork(stats_to_include)",
      details =
        "All node-level network properties are now automatically computed.",
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(cluster_stats)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "buildPublicClusterNetwork(cluster_stats)",
      details =
        "All cluster-level network properties are now automatically computed.",
      env = env,
      user_env = user_env
    )
  }
}

.checkDeprecated.findAssociatedSeqs <- function(
    sample_ids, groups,
    env = caller_env(),
    user_env = caller_env(2)
) {
  if (lifecycle::is_present(groups)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "findAssociatedSeqs(groups)",
      details =
        "Group labels are determined from the unique group ID values.",
      env = env,
      user_env = user_env
    )
  }
  if (lifecycle::is_present(sample_ids)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "findAssociatedSeqs(sample_ids)",
      details = "Custom sample IDs are not relevant to this function.",
      env = env,
      user_env = user_env
    )
  }
}
.checkDeprecated.findAssociatedSeqs2 <- function(
    groups,
    env = caller_env(),
    user_env = caller_env(2)
) {
  if (lifecycle::is_present(groups)) {
    lifecycle::deprecate_warn(
      when = "0.0.9038",
      what = "findAssociatedSeqs(groups)",
      details =
        "Group labels are determined from the unique group ID values.",
      env = env,
      user_env = user_env
    )
  }
}