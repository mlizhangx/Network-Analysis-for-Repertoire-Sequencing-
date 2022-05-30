library(RepSeqNetworkAnalysis)

# load data
dir_main <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
dir_data <- file.path(dir_main, "data")
data <- read.table(
  file.path(dir_data, "TRB-Pt-10-2-500ng-15-04-2020-gDNA_S48.clones.txt"),
  sep = '\t', header = TRUE)



# Default Inputs ----------------------------------------------------------

data = data      # data frame containing req-seq data
nucleo_col = "nSeqCDR3"
amino_col = "aaSeqCDR3"
count_col = "cloneCount"
freq_col = "cloneFraction"
vgene_col = "bestVGene"
dgene_col = "bestDGene"
jgene_col = "bestJGene"
cdr3length_col = "lengthOfCDR3"
other_cols = NULL # other cols to keep; ignored if aggregate_reads = TRUE

# Clone Sequence Settings
clone_seq_type = "amino acid" # or "nucleotide"
min_seq_length = 3 # min clone seq length
drop_chars = NULL # regular expression, e.g. "[*|_]"
aggregate_reads = TRUE
grouping_cols = NULL

# Network Settings
dist_type = "hamming" # or "levenshtein", "hamming", "euclidean_on_atchley"
edge_dist = 1 # max dist for edges

# Network Stats
node_stats = FALSE
stats_to_include = node_stat_settings()
cluster_stats = FALSE

# Plot Settings
plot_title = paste("RepSeq network on CDR3", clone_seq_type, "sequence")
plot_subtitle = ifelse(dist_type == "euclidean_on_atchley",
                       yes = paste("Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of", edge_dist, "between embedded values\n"),
                       no = paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n"))
# size_nodes_by = count_col # can use a double, e.g., 1.0, for fixed size
size_nodes_by = "cloneFraction" # can use a double, e.g., 1.0, for fixed size
node_size_limits = NULL # numeric, length 2
custom_size_legend = NULL
color_nodes_by = "cloneFraction" # use NULL to automatically determine
# color_nodes_by = NULL # use NULL to automatically determine
color_scheme = "default"
custom_color_legend = NULL

# Output Settings
output_dir = file.path(dir_main, "buildRepSeqNetwork", "output") # if NULL, output is not saved to file
save_all = FALSE # by default, only save pdf of plot and csv of node data
data_outfile = "node_data.csv"
plot_outfile = "network_graph.pdf"
plot_width = 12 # passed to pdf()
plot_height = 10 # passed to pdf()
cluster_outfile = "cluster_info.csv" # only saved if save_all = TRUE
igraph_outfile = "network_edgelist.txt" # .txt
matrix_outfile = ifelse(dist_type == "euclidean_on_atchley",
                        yes = "adjacency_matrix.csv",
                        no = "adjacency_matrix.mtx") # .mtx (.csv for euclidean on atchley)
return_all = FALSE # if false, only the node data is returned (unless cluster_stats = TRUE and output_dir = NULL, in which case a list containing the node data and cluster info is returned, with a warning)




# Function ----------------------------------------------------------------


#### INPUT CHECKS ####

# Atchley factor embedding only applicable to amino acid sequences
if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
  stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

# aggregate_reads only applicable to amino acid sequences

# each elem of color_nodes_by is a data column or a node stat to be added


#### PREPARE WORKING ENVIRONMENT ####
# Create output directory if applicable
if (!is.null(output_dir)) { .createOutputDir(output_dir) }

# return type
return_type <- ifelse(return_all,
                      yes = "all",
                      no = ifelse(cluster_stats & is.null(output_dir),
                                  yes = "node_and_cluster_data",
                                  no = "node_data_only"))

# Convert input columns to character if not already
if (is.numeric(nucleo_col)) { nucleo_col <- names(data)[nucleo_col] }
if (is.numeric(amino_col)) { amino_col <- names(data)[amino_col] }
if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }
if (!is.null(color_nodes_by)) {
  if (is.numeric(color_nodes_by)) {
    color_nodes_by <- names(data)[color_nodes_by] } }
if (is.integer(size_nodes_by)) { size_nodes_by <- names(data)[size_nodes_by] }
if (!aggregate_reads) {
  if (is.numeric(vgene_col)) { vgene_col <- names(data)[vgene_col] }
  if (is.numeric(dgene_col)) { dgene_col <- names(data)[dgene_col] }
  if (is.numeric(jgene_col)) { jgene_col <- names(data)[jgene_col] }
  if (is.numeric(other_cols)) { other_cols <- names(data)[other_cols] }
} else {
  if (is.numeric(grouping_cols)) {
    grouping_cols <- names(data)[grouping_cols] }
}

# Add columns in color_nodes_by to other_cols if not present
other_cols <- unique(c(other_cols, color_nodes_by))

# Designate amino acid or nucleotide for clone sequence
clone_seq_col <- amino_col
if (clone_seq_type == "nucleotide") { clone_seq_col <- nucleo_col }


#### FORMAT AND FILTER DATA ####
# Remove sequences below specified length
data <- filterDataBySequenceLength(data, clone_seq_col,
                                   min_length = min_seq_length)
if (nrow(data) < 2) { stop("fewer than two clone sequences meet the specified minimum length") }

# Drop sequences with specified chars
if (!is.null(drop_chars)) {
  data <- data[-grep(drop_chars, data[ , clone_seq_col]), ] }

if (aggregate_reads) { # Aggregate the counts if specified
  data <- aggregateReads(data, clone_seq_col,
                         count_col, freq_col, grouping_cols)

  # Update column name references for count and freq
  if (size_nodes_by == count_col) { size_nodes_by <- "aggCloneCount" }
  if (size_nodes_by == freq_col) { size_nodes_by <- "aggCloneFreq" }
  if (count_col %in% color_nodes_by) {
    color_nodes_by[which(color_nodes_by == count_col)] <- "aggCloneCount" }
  if (freq_col %in% color_nodes_by) {
    color_nodes_by[which(color_nodes_by == freq_col)] <- "aggCloneFreq" }
  count_col <- "aggCloneCount"
  freq_col <- "aggCloneFreq"

} else { # Copy the relevant columns from the input data
  data <-
    data[ , # Keep only the relevant columns:
          c(nucleo_col, amino_col, count_col, freq_col,
            vgene_col, dgene_col, jgene_col, cdr3length_col, other_cols)] }


#### BUILD NETWORK ####
# Generate adjacency matrix for network
adjacency_matrix <-
  generateNetworkFromClones(data[ , clone_seq_col],
                            dist_type, edge_dist,
                            contig_ids = rownames(data),
                            return_type = "adjacency_matrix")

# Subset data to keep only those clones in the network (nonzero degree)
if (dist_type != "euclidean_on_atchley") {
  data <- data[dimnames(adjacency_matrix)[[1]], ] }

# Generate network from adjacency matrix
net <- generateNetworkFromAdjacencyMat(adjacency_matrix)


#### NODE/CLUSTER STATS ####
# Add node-level network characteristics
if (node_stats) {
  data <- addNodeNetworkStats(data, net, stats_to_include) }

# Compute cluster-level network characteristics
if (cluster_stats) {

  if (!"cluster_id" %in% names(data)) {
    data <- addClusterMembership(data, net) }

  degree_col <- NULL
  if ("degree" %in% names(data)) { degree_col <- "degree" }

  cluster_info <- getClusterStats(data, adjacency_matrix, clone_seq_col,
                                  count_col, "cluster_id", degree_col) }


### PLOT(S) OF NETWORK GRAPH ####
# if color_nodes_by is NULL, determine default color variable
if (is.null(color_nodes_by)) {
  if ("degree" %in% names(data)) { # use network degree if available
    color_nodes_by <- "degree"
  } else { # if degree unavailable, color the nodes by clone count
    color_nodes_by <- count_col } }

# Add cluster ID to node color variables if applicable
# if ("cluster_id" %in% names(data) & !"cluster_id" %in% color_nodes_by) {
#   color_nodes_by <- c(color_nodes_by, "cluster_id") }

# Vector of legend titles as node variable names
color_legend_title <- color_nodes_by

# Force consistent legend labels for count/freq if used for node colors
# (ensures size and color share same legend if using same variable)
# if (count_col %in% color_nodes_by) {
#   if (aggregate_reads) {
#     color_legend_title[which(color_legend_title == count_col)] <-
#       "agg clone count"
#   } else {
#     color_legend_title[which(color_legend_title == count_col)] <-
#       "clone count" } }
# if (freq_col %in% color_nodes_by) {
#   if (aggregate_reads) {
#     color_legend_title[which(color_legend_title == count_col)] <-
#       "agg clone freq" }  }

# If multiple coloring variables, extend color scheme to vector if needed
if (length(color_nodes_by) > 1 & length(color_scheme) == 1) {
  color_scheme <- rep(color_scheme, length(color_nodes_by)) }

# Ensure size_nodes_by is a vector or fixed value to use for node sizes
# also ensure size legend title matches color legend title for same variable
size_legend_title <- NULL # default for fixed node size
if (is.character(size_nodes_by)) {
  # if (aggregate_reads) {
  #   if (size_nodes_by == "aggCloneCount") {
  #     size_legend_title <- "agg clone count"
  #   } else if (size_nodes_by == "aggCloneFreq") {
  #     size_legend_title <- "agg clone freq" }
  # } else if (size_nodes_by == count_col) {
  #   size_legend_title <- "clone count"
  # } else {
  size_legend_title <- size_nodes_by
  # }

  size_nodes_by <- data[ , size_nodes_by]
}

# Override legend titles with custom values if provided
if (!is.null(custom_color_legend)) {
  color_legend_title <- custom_color_legend }
if (!is.null(custom_size_legend)) {
  size_legend_title <- custom_size_legend }

# Create one plot for each variable used to color the nodes
temp_plotlist <- list()
for (j in 1:length(color_nodes_by)) {
  cat(paste0("Generating graph plot with nodes colored by ",
             color_nodes_by[[j]], "..."))
  temp_plotlist$newplot <-
    plotNetworkGraph(
      net, title = plot_title, subtitle = plot_subtitle,
      color_nodes_by = data[ , color_nodes_by[[j]]],
      size_nodes_by = size_nodes_by,
      color_legend_title = color_legend_title[[j]],
      size_legend_title = size_legend_title,
      color_scheme = color_scheme[[j]],
      node_size_limits = node_size_limits)
  print(temp_plotlist$newplot) # print to R
  names(temp_plotlist)[[length(names(temp_plotlist))]] <- color_nodes_by[[j]]
  cat(" Done.\n") }


#### SAVE RESULTS ####
if (!aggregate_reads) {  # Rename data columns
  colnames(data)[1:8] <- c(
    "nucleotideSeq", "aminoAcidSeq", "cloneCount", "cloneFrequency",
    "VGene", "DGene", "JGene", "CDR3Length") }

# Save node [& cluster] data
if (!is.null(output_dir)) {
  if (!is.null(data_outfile)) {
    utils::write.csv(data, file = file.path(output_dir, data_outfile),
                     row.names = FALSE)
    cat(paste0("Node-level data saved to file:\n  ", data_outfile, "\n")) }
  if (cluster_stats) {
    utils::write.csv(cluster_info, file = file.path(output_dir,
                                                    cluster_outfile))
    cat(paste0("Cluster-level data saved to file:\n  ",
               file.path(output_dir, cluster_outfile), "\n")) } }

# Save plots to a single pdf
if (!is.null(output_dir) & !is.null(plot_outfile)) {
  grDevices::pdf(file = file.path(output_dir, plot_outfile),
                 width = plot_width, height = plot_height)
  for (j in 1:length(color_nodes_by)) { print(temp_plotlist[[j]]) }
  grDevices::dev.off()
  cat(paste0("Network graph plot saved to file:\n  ",
             file.path(output_dir, plot_outfile), "\n")) }
# Save igraph
if (!is.null(output_dir) & save_all & !is.null(igraph_outfile)) {
  igraph::write_graph(net,
                      file = file.path(output_dir, igraph_outfile),
                      format = "edgelist")
  cat(paste0("Network igraph saved in edgelist format to file:\n  ",
             file.path(output_dir, igraph_outfile), "\n")) }

# Save adjacency matrix
if (!is.null(output_dir) & save_all & !is.null(matrix_outfile)) {
  if (dist_type == "euclidean_on_atchley") {
    utils::write.csv(adjacency_matrix,
                     file.path(output_dir, matrix_outfile),
                     row.names = FALSE)
    cat(paste0("Adjacency matrix saved to file:\n  ",
               file.path(output_dir, matrix_outfile), "\n"))
  } else {
    Matrix::writeMM(adjacency_matrix,
                    file.path(output_dir, matrix_outfile))
    cat(paste0("Adjacency matrix saved to file:\n  ",
               file.path(output_dir, matrix_outfile), "\n")) } }


#### RETURN OUTPUT ####
if (return_type == "node_data_only") {

  cat("All tasks complete.\n")
  # return(data)

} else {
  out <- list("node_data" = data)
  if (cluster_stats) { out$cluster_info <- cluster_info }
  if (return_type == "all") {
    if (length(temp_plotlist == 1)) { out$plot <- temp_plotlist[[1]]
    } else { out$plot <- temp_plotlist[[1]] }
    out$adjacency_matrix <- adjacency_matrix
    out$igraph <- net }
  cat("All tasks complete.\n")
  # return(out)
}
