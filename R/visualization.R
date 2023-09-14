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


addPlots <- function(
    net,
    print_plots = FALSE,
    plot_title = NULL,
    plot_subtitle = "auto",
    color_nodes_by = NULL,
    color_scheme = "default",
    color_legend = "auto",
    color_title = "auto",
    edge_width = 0.1,
    size_nodes_by = 0.5,
    node_size_limits = NULL,
    size_title = "auto",
    verbose = FALSE
) {
  .MUST.isBaseNetworkOutput(net)
  if (.hasPlots(net)) {
    .MUST.isNamedList(net$plots)
    .checkPlotsAgainstLayout(net$plots)
  }
  print_plots <- .checkTF(print_plots, TRUE)
  plot_title <- .check(plot_title, .isString, NULL, ornull = TRUE)
  plot_subtitle <- .check(plot_subtitle, .isString, "auto", ornull = TRUE)
  color_nodes_by <- .checkColorNodesBy(color_nodes_by, net$node_data)
  color_scheme <- .checkColorScheme(color_scheme, color_nodes_by, "default")
  color_legend <- .checkTF(color_legend, "auto", or_auto = TRUE)
  color_title <- .checkColorTitle(color_title, color_nodes_by)
  edge_width <- .check(edge_width, .isPos, 0.1)
  size_nodes_by <- .checkSizeNodesBy(size_nodes_by, net$node_data)
  node_size_limits <- .checkNodeSizeLimits(node_size_limits)
  size_title <- .check(size_title, .isString, "auto", ornull = TRUE)
  msg <- .makemsg(verbose)
  if (isTRUE("details" %in% names(net))) {
    plot_subtitle <- .makePlotSubtitle(
      plot_subtitle,
      seq_col = net$details$seq_col,
      dist_type = net$details$dist_type,
      dist_cutoff = net$details$dist_cutoff
    )
  } else if (isTRUE(plot_subtitle == "auto")) {
    plot_subtitle <- NULL
  }
  if (isTRUE("plots" %in% names(net))) {
    if ("graph_layout" %in% names(net$plots)) {
      layout <- net$plots$graph_layout
    } else {
      layout <- as.matrix(net$plots[[1]]$data[c("x", "y")])
    }
    newplots <- .generateNetworkGraphPlots(
      net$igraph, net$node_data,
      print_plots, plot_title, plot_subtitle, color_nodes_by,
      color_scheme, color_legend, color_title, edge_width, size_nodes_by,
      node_size_limits, size_title, layout, verbose
    )
    # avoid duplicate element names:
    net$plots <- net$plots[!names(net$plots) %in% names(newplots)]
    net$plots <- c(net$plots, newplots)
    msg("New plots added to ", sQuote("net$plots"))
  } else {
    net$plots <- .generateNetworkGraphPlots(
      net$igraph, net$node_data,
      print_plots, plot_title, plot_subtitle, color_nodes_by,
      color_scheme, color_legend, color_title, edge_width, size_nodes_by,
      node_size_limits, size_title, verbose = verbose
    )
    msg("New plots added to ", sQuote("net$plots"))
  }
  net
}


generateNetworkGraphPlots <- function(
    igraph,
    data,
    print_plots = FALSE,
    plot_title = NULL,
    plot_subtitle = NULL,
    color_nodes_by = NULL,
    color_scheme = "default",
    color_legend = "auto",
    color_title = "auto",
    edge_width = 0.1,
    size_nodes_by = 0.5,
    node_size_limits = NULL,
    size_title = "auto",
    layout = NULL,
    verbose = FALSE
) {
  .MUST.isIgraph(igraph)
  data_name <- deparse(substitute(data))
  data <- as.data.frame(data)
  .MUST.isDataFrame(data, data_name)
  .MUST.doesIgraphMatchData(igraph, data)
  layout <- .check(layout, .doesLayoutMatchData, NULL, ornull = TRUE,
                   data = data
  )
  print_plots <- .checkTF(print_plots, TRUE)
  plot_title <- .check(plot_title, .isString, NULL, ornull = TRUE)
  plot_subtitle <- .check(plot_subtitle, .isString, "auto", ornull = TRUE)
  color_nodes_by <- .checkColorNodesBy(color_nodes_by, data)
  color_scheme <- .checkColorScheme(color_scheme, color_nodes_by, "default")
  color_legend <- .checkTF(color_legend, "auto", or_auto = TRUE)
  color_title <- .checkColorTitle(color_title, color_nodes_by)
  edge_width <- .check(edge_width, .isPos, 0.1)
  size_nodes_by <- .checkSizeNodesBy(size_nodes_by, data)
  node_size_limits <- .checkNodeSizeLimits(node_size_limits)
  size_title <- .check(size_title, .isString, "auto", ornull = TRUE)
  .generateNetworkGraphPlots(
    igraph, data,
    print_plots, plot_title, plot_subtitle, color_nodes_by,
    color_scheme, color_legend, color_title, edge_width, size_nodes_by,
    node_size_limits, size_title, layout, verbose
  )
}


labelNodes <- function(net, node_labels, plots = NULL,
                       size = 5, color = "black"
) {
  .MUST.isBaseNetworkOutput(net, "net")
  .MUST.hasPlots(net, "net")
  if (is.null(plots)) {
    plots <- names(net$plots)[names(net$plots) != "graph_layout"]
  }
  .MUST.isCharOrPosIntegerVector(plots, "plots")
  node_labels <- as.vector(node_labels)
  .MUST.isCharOrNumericVector(node_labels, "node_labels")
  if (length(node_labels) == 1) {
    .stopifnot(.isDataColref(node_labels, net$node_data),
               "node_labels", "has length 1",
               "but is not a valid column reference for 'net$node_data'"
    )
    node_labels <- net$node_data[[node_labels]]
  }
  size <- .check(size, .isPos, 5)
  color <- .check(color, .isString, "black")
  for (i in 1:length(plots)) {
    if (
      (is.character(plots) && .hasElem(net$plots, plots[[i]]) ||
       is.numeric(plots) && plots[[i]] <= length(net$plots)
      ) &&
      .isGgraph(net$plots[[plots[[i]]]]) &&
      .hasLength(node_labels, length(net$plots[[plots[[i]]]]$data$x))
    ) {
      net$plots[[plots[[i]]]] <- net$plots[[plots[[i]]]] +
        ggraph::geom_node_text(
          ggplot2::aes(label = node_labels),
          size = size, color = color
        )
    } else {
      arg <- plots[[i]]
      if (is.character(arg)) { arg <- dQuote(arg) }
      warning("entry ", i, " of 'plots' (", arg, ") does not reference a ",
              "valid ggraph in 'net$plots'"
      )
    }
  }
  net
}


addGraphLabels <- function(plot, node_labels, size = 5, color = "black") {
  .MUST.isGgraph(plot, "plot")
  node_labels <- as.vector(node_labels)
  .MUST.isCharOrNumericVector(node_labels)
  .MUST.hasLength(node_labels, length(plot$data$x))
  size <- .check(size, .isPos, 5)
  color <- .check(color, .isString, "black")
  plot +
    ggraph::geom_node_text(
      ggplot2::aes(label = node_labels),
      size = size, color = color
    )
}

labelClusters <- function(
    net,
    plots = NULL,
    top_n_clusters = 20,
    cluster_id_col = "cluster_id",
    criterion = "node_count",
    size = 5, color = "black",
    greatest_values = TRUE
) {
  .MUST.isBaseNetworkOutput(net, "net")
  .MUST.hasPlots(net, "net")
  if (is.null(plots)) {
    plots <- names(net$plots)[names(net$plots) != "graph_layout"]
  }
  .MUST.isCharOrPosIntegerVector(plots, "plots")
  .MUST.isString(cluster_id_col)
  .MUST.isDataColref(cluster_id_col, net$node_data)
  .MUST.isIntegerVector(net$node_data[[cluster_id_col]], factor_ok = TRUE)
  top_n_clusters <- .check(top_n_clusters, .isPosInt, 20)
  size <- .check(size, .isPos, 5)
  color <- .check(color, .isString, "black")
  greatest_values <- .checkTF(greatest_values, TRUE)
  criterion <- .check(criterion, .isString, "node_count")
  has_cdat <- .hasClusterData(net)
  cdat_matches <- has_cdat && .hasDetails(net) &&
    .hasElem(net$details, "cluster_data_goes_with") &&
    net$details$cluster_data_goes_with == cluster_id_col
  if (criterion != "node_count" &&
      (!has_cdat || !cdat_matches ||
       !.isDataColref(criterion, net$cluster_data)
      )
  ) {
    warning(dQuote(criterion), " invalid for ", sQuote("criterion"),
            " since ", sQuote("net$cluster_data"),
            " does not contain a column with this name ",
            "or does not correspond to cluster membership variable ",
            dQuote(cluster_id_col), " in ", sQuote("net$node_data"), ". ",
            "Defaulting to ", dQuote("node_count")
    )
    criterion <- "node_count"
  }
  if (cdat_matches) {
    cdat <- net$cluster_data
  } else { # compute node count for ranking
    cdat <- as.data.frame(table(net$node_data[[cluster_id_col]]))
    names(cdat) <- c("cluster_id", "node_count")
    criterion <- "node_count"
  }
  if (isTRUE(greatest_values)) {
    cdat <- cdat[order(-cdat[[criterion]]) , ]
  } else {
    cdat <- cdat[order(cdat[[criterion]]) , ]
  }
  clusters <- as.integer(cdat$cluster_id[1:top_n_clusters])
  node_labels <- as.integer(net$node_data[[cluster_id_col]])
  node_labels[duplicated(node_labels) | !node_labels %in% clusters] <- NA
  labelNodes(net, node_labels, plots, size, color)
}


addClusterLabels <- function(
    plot,
    net,
    top_n_clusters = 20,
    cluster_id_col = "cluster_id",
    criterion = "node_count",
    size = 5, color = "black",
    greatest_values = TRUE
) {
  .MUST.isGgraph(plot, "plot")
  .MUST.isBaseNetworkOutput(net, "net")
  .MUST.doesPlotMatchData(plot, net$node_data)
  .MUST.isString(cluster_id_col)
  .MUST.isDataColref(cluster_id_col, net$node_data)
  .MUST.isIntegerVector(net$node_data[[cluster_id_col]], factor_ok = TRUE)
  top_n_clusters <- .check(top_n_clusters, .isPosInt, 20)
  size <- .check(size, .isPos, 5)
  color <- .check(color, .isString, "black")
  greatest_values <- .checkTF(greatest_values, TRUE)
  criterion <- .check(criterion, .isString, "node_count")
  has_cdat <- .hasClusterData(net)
  cdat_matches <- has_cdat && .hasDetails(net) &&
    .hasElem(net$details, "cluster_data_goes_with") &&
    net$details$cluster_data_goes_with == cluster_id_col
  if (criterion != "node_count" &&
      (!has_cdat || !cdat_matches ||
       !.isDataColref(criterion, net$cluster_data)
      )
  ) {
    warning(dQuote(criterion), " invalid for ", sQuote("criterion"),
            " since ", sQuote("net$cluster_data"),
            " does not contain a column with this name ",
            "or does not correspond to cluster membership variable ",
            dQuote(cluster_id_col), " in ", sQuote("net$node_data"), ". ",
            "Defaulting to ", dQuote("node_count")
    )
    criterion <- "node_count"
  }
  if (cdat_matches) {
    cdat <- net$cluster_data
  } else { # compute node count for ranking
    cdat <- as.data.frame(table(net$node_data[[cluster_id_col]]))
    names(cdat) <- c("cluster_id", "node_count")
    criterion <- "node_count"
  }
  if (isTRUE(greatest_values)) {
    cdat <- cdat[order(-cdat[[criterion]]) , ]
  } else {
    cdat <- cdat[order(cdat[[criterion]]) , ]
  }
  clusters <- as.integer(cdat$cluster_id[1:top_n_clusters])
  node_labels <- as.integer(net$node_data[[cluster_id_col]])
  node_labels[duplicated(node_labels) | !node_labels %in% clusters] <- NA
  addGraphLabels(plot, node_labels, size, color)
}


extractLayout <- function(plot) {
  .MUST.isGgraph(plot)
  as.matrix(plot$data[c("x", "y")])
}


plotNetworkGraph <- function(
    igraph,
    plot_title = NULL,
    plot_subtitle = NULL,
    color_nodes_by = NULL,
    color_scheme = "default",
    color_legend = "auto",
    color_title = "auto",
    edge_width = 0.1,
    size_nodes_by = 0.5,
    node_size_limits = NULL,
    size_title = "auto",
    outfile = NULL,
    pdf_width = 12,
    pdf_height = 8
) {
  lifecycle::deprecate_soft(
    when = "1.0.1",
    what = "plotNetworkGraph()",
    with = "addPlots()"
  )
  graph_plot <-
    .plotNetworkGraph(
      igraph, NULL,
      plot_title,
      plot_subtitle,
      color_nodes_by,
      color_scheme,
      color_legend,
      color_title,
      edge_width,
      size_nodes_by,
      node_size_limits,
      size_title
    )
  if (!is.null(outfile)) {
    grDevices::pdf(file = outfile, width = pdf_width, height = pdf_height)
    print(graph_plot)
    grDevices::dev.off()
  }
  graph_plot
}

# Internal ----------------------------------------------------------------

.generateNetworkGraphPlots <- function(
    igraph,
    data,
    print_plots = TRUE,
    plot_title = NULL,
    plot_subtitle = NULL,
    color_nodes_by = NULL,
    color_scheme = "default",
    color_legend = "auto",
    color_title = "auto",
    edge_width = 0.1,
    size_nodes_by = 0.5,
    node_size_limits = NULL,
    size_title = "auto",
    layout = NULL,
    verbose = FALSE
) {
  msg <- .makemsg(verbose)
  if (is.character(size_nodes_by) && !is.numeric(data[[size_nodes_by]])) {
    warning("Non-numeric variable specified for ", sQuote("size_nodes_by"),
            ". Defaulting to fixed node size"
    )
    size_nodes_by <- 0.5
  }
  new_inputs <- .harmonizePlottingInputs(
    color_nodes_by, color_scheme, color_title, size_nodes_by, size_title
  )
  color_scheme <- new_inputs$color_scheme
  color_title <- new_inputs$color_title
  size_title <- new_inputs$size_title
  if (is.character(size_nodes_by)) {
    size_nodes_by <- data[[size_nodes_by]]
  }
  plotlist <- list()
  if (is.null(layout)) { layout <- igraph::layout_components(igraph) }
  if (is.null(color_nodes_by)) {
    msg("Generating graph plot...", newline = FALSE)
    plotlist$uniform_color <- .plotNetworkGraph(
      igraph, layout,
      plot_title = plot_title,
      plot_subtitle = plot_subtitle,
      color_nodes_by = NULL,
      color_title = NULL,
      color_scheme = color_scheme,
      color_legend = FALSE,
      edge_width = edge_width,
      size_nodes_by = size_nodes_by,
      size_title = size_title,
      node_size_limits = node_size_limits
    )
    if (print_plots) { print(plotlist$uniform_color) }
    msg(" Done.")
  } else {
    for (j in 1:length(color_nodes_by)) {
      msg("Generating graph plot with nodes colored by ", color_nodes_by[[j]],
          "...", newline = FALSE
      )
      plotlist$newplot <- .plotNetworkGraph(
        igraph, layout,
        plot_title = plot_title,
        plot_subtitle = plot_subtitle,
        color_nodes_by = data[[color_nodes_by[[j]]]],
        color_title = color_title[[j]],
        color_scheme = color_scheme[[j]],
        color_legend = color_legend,
        edge_width = edge_width,
        size_nodes_by = size_nodes_by,
        size_title = size_title,
        node_size_limits = node_size_limits
      )
      if (print_plots) { print(plotlist$newplot) }
      names(plotlist)[[length(names(plotlist))]] <- color_nodes_by[[j]]
      msg(" Done.")
    }
  }
  plotlist$graph_layout <- layout
  plotlist
}


.plotNetworkGraph <- function(
    igraph, layout = NULL,
    plot_title = NULL,
    plot_subtitle = NULL,
    color_nodes_by = NULL,
    color_scheme = "default",
    color_legend = "auto",
    color_title = "auto",
    edge_width = 0.1,
    size_nodes_by = 0.5,
    node_size_limits = NULL,
    size_title = "auto"
) {

  # Default legend titles
  if (length(size_nodes_by) > 1 && isTRUE(size_title == "auto")) {
    size_title <- deparse(substitute(size_nodes_by))
  }
  if (!is.null(color_nodes_by) && isTRUE(color_title == "auto")
  ) {
    color_title <- deparse(substitute(color_nodes_by))
  }

  # Is color variable continuous or discrete?
  color_type <- "continuous"
  if (!is.null(color_nodes_by)) {
    color_type <- ggplot2::scale_type(color_nodes_by)[[1]]
  }

  # Is the color legend too big to show?
  if (isTRUE(color_legend == "auto")) {
    color_legend <- FALSE
    too_many_legend_values <-
      color_type == "discrete" && length(unique(color_nodes_by)) > 20
    if (!is.null(color_nodes_by) && !too_many_legend_values) {
      color_legend <- TRUE
    }
    if (too_many_legend_values) {
      if (!is.null(plot_subtitle)) {
        plot_subtitle <- paste0(plot_subtitle, "\n")
      }
      plot_subtitle <- paste0(plot_subtitle, "Nodes colored by ", color_title)
    }
  }

  # Create base plot
  if (is.null(layout)) { layout <- igraph::layout_components(igraph) }
  graph_plot <-
    ggraph::ggraph(igraph, layout = layout) +
    ggraph::geom_edge_link0(width = edge_width, colour = "grey")

  graph_plot <- .addGeomNodePoint(
    graph_plot, color_nodes_by, size_nodes_by, node_size_limits
  ) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(title = plot_title,
                  subtitle = plot_subtitle
    )

  # Legend and color scale
  if (length(size_nodes_by) > 1) {
    graph_plot <- .setSizeLegendTitle(graph_plot, size_title)
  }
  if (!is.null(color_nodes_by)) {
    graph_plot <- .setColorLegendTitle(graph_plot, color_legend, color_title,
                                       color_type, size_nodes_by,
                                       color_nodes_by, size_title
    )
    graph_plot <- .addColorScheme(graph_plot, color_type, color_scheme,
                                  n_colors = length(unique(color_nodes_by))
    )
  }
  graph_plot
}



.harmonizePlottingInputs <- function(
    color_nodes_by, color_scheme, color_title, size_nodes_by, size_title)
{
  if (is.null(color_nodes_by)) {
    color_title <- NULL
  } else { # color_nodes_by is non NULL
    if (length(color_nodes_by) > 1) { # check for duplicate entries
      if (sum(duplicated(color_nodes_by)) > 0) {
        duplicates <- duplicated(color_nodes_by)
        color_nodes_by <- color_nodes_by[!duplicates]
        if (length(color_scheme) > 1) {
          color_scheme <- color_scheme[!duplicates]
          warning(
            "duplicate entries removed from ", sQuote("color_nodes_by"),
            ";\n corresponding entries removed from ", sQuote("color_scheme")
          )
        } else {

          warning("Duplicate entries removed from ", sQuote("color_nodes_by"))
        }
      }
    }
    if (length(color_nodes_by) > 1) { # might need to extend scalars to vectors
      if (length(color_scheme) == 1) {
        color_scheme <- rep(color_scheme, length(color_nodes_by))
      }
      if (!is.null(color_title)) {
        if (length(color_title) == 1) {
          color_title <- rep(color_title, length(color_nodes_by))
        }
      } else { # color_title is NULL
        color_title <- rep("", length(color_nodes_by))
      }
    }
    if (!is.null(color_title)) {
      for (i in 1:length(color_title)) {
        if (isTRUE(color_title[[i]] == "auto")) {
          color_title[[i]] <- color_nodes_by[[i]]
        }
      }
    }
  }
  if (!is.null(size_title) && isTRUE(size_title == "auto")) {
    if (is.numeric(size_nodes_by)) { size_title <- NULL }
    if (is.character(size_nodes_by)) { size_title <- size_nodes_by }
  }
  list(color_scheme = color_scheme,
       color_title = color_title,
       size_title = size_title
  )
}

.passColorNodesBy <- function(color_nodes_by, data, count_col) {

  if (length(color_nodes_by) == 1 && isTRUE(color_nodes_by == "auto")) {
    if (!is.null(count_col)) {
      color_nodes_by <- count_col
    } else if ("cluster_id" %in% names(data)) {
      color_nodes_by <- "cluster_id"
    } else if ("degree" %in% names(data)) {
      color_nodes_by <- "degree"
    } else if ("transitivity" %in% names(data)) {
      color_nodes_by <- "transitivity"
    } else if ("closeness" %in% names(data)) {
      color_nodes_by <- "closeness"
    } else if ("eigen_centrality" %in% names(data)) {
      color_nodes_by <- "eigen_centrality"
    } else if ("betweenness" %in% names(data)) {
      color_nodes_by <- "betweenness"
    } else if ("coreness" %in% names(data)) {
      color_nodes_by <- "coreness"
    } else if ("authority_score" %in% names(data)) {
      color_nodes_by <- "authority_score"
    } else if ("page_rank" %in% names(data)) {
      color_nodes_by <- "page_rank"
    } else {
      color_nodes_by <- NULL
    }
  }
  color_nodes_by

}

.makePlotTitle <- function(plot_title, type = "standard", network_name = NULL) {

  if (is.null(plot_title)) {
    return(plot_title)
  }
  if (type == "standard") {
    if (isTRUE(plot_title == "auto")) {
      if (!is.null(network_name)) {
        return(network_name)
      } else {
        return("Network by Receptor Sequence Similarity")
      }
    }
  } else if (type == "pub_clust_rep") {
    if (isTRUE(plot_title == "auto")) {
      return("Public Clusters Using Representative TCR/BCRs")
    }
  }
  plot_title

}

.makePlotSubtitle <- function(plot_subtitle, type = "standard", seq_col = NULL,
                              dist_type = NULL, dist_cutoff = NULL
) {
  if (is.null(plot_subtitle)) {
    return(plot_subtitle)
  }
  if (type == "standard") {
    if (isTRUE(plot_subtitle == "auto")) {
      if (length(seq_col) == 2) {
        return(paste(
          "Each node denotes a single TCR/BCR cell",
          "\nEdges denote a maximum", dist_type, "distance of", dist_cutoff,
          "between sequences in corresponding chains\n"
        ))
      }
      return(paste(
        "Each node denotes a single TCR/BCR cell or clone",
        "\nEdges denote a maximum", dist_type, "distance of", dist_cutoff,
        "between receptor sequences\n"
      ))
    }
  } else if (type == "pub_clust_rep") {
    if (plot_subtitle == "auto" && seq_col == "seq_w_max_count") {
      return(paste("Each node represents the most abundant TCR/BCR sequence",
                   "\nfrom a distinct cluster within a single sample",
                   "\nEdges denote a maximum", dist_type, "distance of",
                   dist_cutoff, "between sequences\n"
      ))
    }
    return(paste("Each node represents a representative TCR/BCR sequence",
                 "\nfrom a distinct cluster within a single sample",
                 "\nRepresentative based on property:", seq_col,
                 "\nEdges denote a maximum", dist_type, "distance of",
                 dist_cutoff, "between sequences\n"
    ))
  }
  plot_subtitle
}


.addGeomNodePoint <- function(
    graph_plot, color_nodes_by, size_nodes_by, node_size_limits
) {

  if (!is.null(color_nodes_by)) {
    if (length(size_nodes_by) > 1) { # colors and sizes are variable
      graph_plot <- graph_plot +
        ggraph::geom_node_point(
          ggplot2::aes(color = color_nodes_by,
                       size = size_nodes_by
          )
        )
    } else if (length(size_nodes_by) == 1) {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(
          ggplot2::aes(color = color_nodes_by),
          size = size_nodes_by
        )
    } else { # size_nodes_by = NULL
      graph_plot <- graph_plot +
        ggraph::geom_node_point(
          ggplot2::aes(color = color_nodes_by)
        )
    }
  } else { # color_nodes_by = NULL
    if (length(size_nodes_by) > 1) {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(
          ggplot2::aes(size = size_nodes_by)
        )
    } else if (length(size_nodes_by) == 1) {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(
          size = size_nodes_by
        )
    } else { # size_nodes_by = NULL
      graph_plot <- graph_plot +
        ggraph::geom_node_point()
    }
  }

  if (length(size_nodes_by) > 1 && length(node_size_limits) == 2) {
    graph_plot <- graph_plot +
      ggplot2::scale_size(range = node_size_limits)
  }

  graph_plot
}


.setSizeLegendTitle <- function(graph_plot, size_title) {
  graph_plot <- graph_plot +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title = size_title
      )
    )
  graph_plot
}


.setColorLegendTitle <- function(graph_plot, color_legend, color_title,
                                 color_type, size_nodes_by, color_nodes_by,
                                 size_title
) {
  if (isTRUE(color_legend)) {
    is_not_combined_legends <-
      !isTRUE(all.equal(size_nodes_by, color_nodes_by)) ||
      isTRUE(color_title != size_title) ||
      isTRUE(sum(c(is.null(color_title), is.null(size_title))) == 1)
    # use colorbar for cts unless color & size use combined legend
    if (color_type == "continuous" && is_not_combined_legends) {
      graph_plot <- graph_plot +
        ggplot2::guides(
          color = ggplot2::guide_colorbar(
            title = color_title
          )
        )
    } else {
      graph_plot <- graph_plot +
        ggplot2::guides(
          color = ggplot2::guide_legend(
            title = color_title
          )
        )
    }
  } else {
    graph_plot <- graph_plot +
      ggplot2::guides(
        color = "none"
      )
  }
  graph_plot
}

.addColorScheme <- function(graph_plot, color_type, color_scheme, n_colors) {

  if (isTRUE(color_scheme == "default")) {
    return(graph_plot)
  }
  viridis_keys <- c(
    "A", "B", "C", "F", "G", "magma", "inferno", "plasma", "rocket", "mako"
  )
  if (color_type == "continuous") {
    if (color_scheme %in% viridis_keys) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = color_scheme, begin = 0.05, end = 0.95
        )
    } else if (color_scheme %in% c("D", "viridis")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = color_scheme, begin = 0, end = 0.95
        )
    } else if (color_scheme %in% c("E", "H", "cividis",  "turbo")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(option = color_scheme)
    } else if (color_scheme %in% paste0(viridis_keys, "-1")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = strsplit(color_scheme, "-1")[[1]],
          begin = 0.05, end = 0.95, direction = -1
        )
    } else if (color_scheme %in% c("D-1", "viridis-1")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = strsplit(color_scheme, "-1")[[1]],
          begin = 0, end = 0.95, direction = -1
        )
    } else if (color_scheme %in% c("E-1", "H-1", "cividis-1",  "turbo-1")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = strsplit(color_scheme, "-1")[[1]], direction = -1
        )
    } else {
      warning("Value for color_scheme is not a valid option for ",
              "continuous variables; using default color scheme instead"
      )
    }
  } else { # discrete color scheme
    if (color_scheme %in% viridis_keys) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = color_scheme, discrete = TRUE, begin = 0.05, end = 0.95
        )
    } else if (color_scheme %in% c("D", "viridis")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = color_scheme, discrete = TRUE, begin = 0, end = 0.95
        )
    } else if (color_scheme %in% c("E", "H", "cividis",  "turbo")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(option = color_scheme, discrete = TRUE)
    } else if (color_scheme %in% paste0(viridis_keys, "-1")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
          begin = 0.05, end = 0.95, direction = -1
        )
    } else if (color_scheme %in% c("D-1", "viridis-1")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = strsplit(color_scheme, "-1")[[1]], discrete = TRUE,
          begin = 0, end = 0.95, direction = -1
        )
    } else if (color_scheme %in% c("E-1", "H-1", "cividis-1",  "turbo-1")) {
      graph_plot <- graph_plot +
        ggraph::scale_color_viridis(
          option = strsplit(color_scheme, "-1")[[1]],
          discrete = TRUE, direction = -1)
    } else if (color_scheme %in% grDevices::hcl.pals()) {
      graph_plot <- graph_plot +
        ggplot2::scale_color_manual(
          values = grDevices::hcl.colors(
            n = n_colors, palette = color_scheme
          )
        )
    } else {
      warning("Value for color_scheme is not a valid option for ",
              "discrete variables; using default color scheme instead"
      )
    }
  }
  graph_plot
}
