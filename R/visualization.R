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
  set.seed(9999)
  layout <- igraph::layout_components(igraph)
  graph_plot <-
    ggraph::ggraph(igraph, layout = layout) +
    ggraph::geom_edge_link0(width = edge_width, colour = "grey")
  if (!is.null(color_nodes_by)) {
    if (!is.null(size_nodes_by)) {
      if (is.numeric(size_nodes_by) && length(size_nodes_by) == 1) {
        graph_plot <- graph_plot +
          ggraph::geom_node_point(
            ggplot2::aes(color = color_nodes_by), size = size_nodes_by
          )
      } else if (length(size_nodes_by) > 1) {
        graph_plot <- graph_plot +
          ggraph::geom_node_point(
            ggplot2::aes(color = color_nodes_by, size = size_nodes_by)
          )
        if (length(node_size_limits) == 2) {
          graph_plot <-
            graph_plot + ggplot2::scale_size(range = node_size_limits)
        }
        if (!is.null(size_title) && size_title == "auto") {
          size_title <- deparse(substitute(size_nodes_by))
        }
      }
    } else { # size_nodes_by is null or invalid
      graph_plot <- graph_plot +
        ggraph::geom_node_point(ggplot2::aes(color = color_nodes_by))

    }
  } else if (!is.null(size_nodes_by)) { # color_nodes_by is null or invalid
    if (is.numeric(size_nodes_by) && length(size_nodes_by) == 1) {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(size = size_nodes_by)

    } else if (length(size_nodes_by) > 1) {
      graph_plot <- graph_plot +
        ggraph::geom_node_point(ggplot2::aes(size = size_nodes_by))
      if (!is.null(node_size_limits)) {
        graph_plot <-
          graph_plot + ggplot2::scale_size(range = node_size_limits)
      }
      if (!is.null(size_title) && size_title == "auto") {
        size_title <- deparse(substitute(size_nodes_by))
      }
    }
  } else { # size_nodes_by null or invalid
    graph_plot <- graph_plot + ggraph::geom_node_point()
  }

  if (is.null(color_nodes_by)) {
    color_type <- "continuous"
  } else {
    color_type <- ggplot2::scale_type(color_nodes_by)[[1]]
  }

  if (color_legend == "auto") {
    if (is.null(color_nodes_by)) {
      color_legend <- FALSE
    } else {
      if (color_type == "continuous" || length(unique(color_nodes_by)) <= 20) {
        color_legend <- TRUE
      } else { # too many legend values
        color_legend <- FALSE
        color_varname <- ifelse(is.null(color_title),
                                no = color_title,
                                yes = deparse(substitute(color_nodes_by))
        )
        subtitle_affix <- paste0("Nodes colored by ", color_varname)
        plot_subtitle <- ifelse(is.null(plot_subtitle),
                                yes = subtitle_affix,
                                no = paste0(plot_subtitle, "\n", subtitle_affix)
        )
      }
    }
  }

  graph_plot <- graph_plot +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::labs(title = plot_title, subtitle = plot_subtitle)

  if (color_legend == TRUE) {
    if (!is.null(color_title) && color_title == "auto") {
      color_title <- deparse(substitute(color_nodes_by))
    }
    graph_plot <- graph_plot +
      ggplot2::guides(color = ggplot2::guide_legend(title = color_title))
  } else {
    graph_plot <- graph_plot + ggplot2::guides(color = "none")
  }
  graph_plot <- graph_plot +
    ggplot2::guides(size = ggplot2::guide_legend(title = size_title))

  if (color_scheme != "default") {
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
        warning(paste("Value for color_scheme is not a valid option for",
                      "continuous variables; using default color scheme instead"

        ))
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
              n = length(unique(color_nodes_by)), palette = color_scheme
            )
          )
      } else {
        warning(paste("Value for color_scheme is not a valid option for",
                      "discrete variables; using default color scheme instead"
        ))
      }
    }
  }
  if (!is.null(outfile)) {
    grDevices::pdf(file = outfile, width = pdf_width, height = pdf_height)
    print(graph_plot)
    grDevices::dev.off()
    cat(paste0("Plot of network graph saved to file:\n  ", outfile, "\n"))
  }
  graph_plot
}


generateNetworkGraphPlots <- function(
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
    size_title = "auto"
) {
  .checkargs.generateNetworkGraphPlots(
    igraph, data, print_plots, plot_title, plot_subtitle, color_nodes_by,
    color_scheme, color_legend, color_title, edge_width, size_nodes_by,
    node_size_limits, size_title
  )
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
  if (is.null(color_nodes_by)) {
    cat("Generating graph plot...")
    plotlist$uniform_color <- plotNetworkGraph(
      igraph,
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
    cat(" Done.\n")
  } else {
    for (j in 1:length(color_nodes_by)) {
      cat(paste0("Generating graph plot with nodes colored by ",
                 color_nodes_by[[j]], "..."))
      plotlist$newplot <- plotNetworkGraph(
        igraph,
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
      cat(" Done.\n")
    }
  }
  plotlist
}

addGraphLabels <- function(plot, node_labels, size = 5, color = "black") {
  plot +
    ggraph::geom_node_text(
      ggplot2::aes(label = node_labels),
      size = size, color = color
    )
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
  .checkargs.addClusterLabels(
    plot, net, top_n_clusters, cluster_id_col, criterion,
    size, color, greatest_values
  )
  dat <- net$node_data
  cdat <- net$cluster_data
  if (isTRUE(greatest_values)) {
    cdat <- cdat[order(-cdat[[criterion]]) , ]
  } else {
    cdat <- cdat[order(cdat[[criterion]]) , ]
  }
  clusters <- as.integer(cdat$cluster_id[1:top_n_clusters])
  node_labels <- as.integer(dat[[cluster_id_col]])
  node_labels[duplicated(node_labels) | !node_labels %in% clusters] <- NA
  addGraphLabels(plot, node_labels, size, color)
}



# Internal ----------------------------------------------------------------
.generateNetworkGraphPlotsGuarded <- function(
    igraph, data, print_plots,
    plot_title = NULL, plot_subtitle = NULL,
    color_nodes_by = NULL, color_scheme = "default",
    color_legend = "auto", color_title = "auto",
    edge_width = 0.1, size_nodes_by = 0.5,
    node_size_limits = NULL, size_title = "auto")
{
  if (nrow(data) > 1e06) {
    warning(paste(
      "Network contains over 1 million nodes; depending on the number of",
      "network edges, this may exceed ggraph limitations. Skipping automatic",
      "generation of network graph; you can attempt to generate the graph",
      "manually using generateNetworkPlots()"
    ))
    return(invisible(NULL))
  } else {
    return(generateNetworkGraphPlots(
      igraph, data, print_plots, plot_title, plot_subtitle, color_nodes_by,
      color_scheme, color_legend, color_title, edge_width, size_nodes_by,
      node_size_limits, size_title
    ))
  }
}

.harmonizePlottingInputs <- function(
    color_nodes_by, color_scheme, color_title, size_nodes_by, size_title)
{
  if (is.null(color_nodes_by)) {
    color_title <- NULL
  } else { # color_nodes_by is non NULL
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
        if (color_title[[i]] == "auto") {
          color_title[[i]] <- color_nodes_by[[i]]
        }
      }
    }
  }
  if (!is.null(size_title) && size_title == "auto") {
    if (is.numeric(size_nodes_by)) { size_title <- NULL }
    if (is.character(size_nodes_by)) { size_title <- size_nodes_by }
  }
  list(color_scheme = color_scheme,
       color_title = color_title,
       size_title = size_title
  )
}

.passColorNodesBy <- function(color_nodes_by, data, count_col) {

  if (length(color_nodes_by) == 1 && color_nodes_by == "auto") {
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
    if (plot_title == "auto") {
      if (!is.null(network_name)) {
        return(network_name)
      } else {
        return("Network by Receptor Sequence Similarity")
      }
    }
  } else if (type == "pub_clust_rep") {
    if (plot_title == "auto") {
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
    if (plot_subtitle == "auto") {
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
      return(paste("Each node represents a TCR/BCR sequence.",
                   "Similar sequences are joined by edges.",
                   "\nNetwork includes one representative TCR/BCR sequence",
                   "from each cluster within each sample",
                   "\nRepresentative sequence from each cluster:",
                   "TCR/BCR with greatest clone count"
      ))
    }
    return(paste("Each node represents a TCR/BCR sequence.",
                 "Similar sequences are joined by edges.",
                 "\nNetwork includes one representative TCR/BCR sequence",
                 "from each cluster within each sample",
                 "\nRepresentative sequence from each cluster",
                 "based on property:", seq_col
    ))
  }
  plot_subtitle
}

