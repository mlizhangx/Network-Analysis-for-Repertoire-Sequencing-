
# Prepare workspace -------------------------------------------------------



library(RepSeqNetworkAnalysis)

# set working directory to dir of current script
dir_main <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
dir_data <- file.path(dir_main, "data")
# load data sets
data <- read.table(
  file.path(dir_data, "TRB-Pt-10-2-500ng-15-04-2020-gDNA_S48.clones.txt"),
  sep = '\t', header = TRUE)


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
cluster_stats = TRUE

# Plot Settings
plot_title = paste("RepSeq network on CDR3", clone_seq_type, "sequence")
plot_subtitle = ifelse(dist_type == "euclidean_on_atchley",
                       yes = paste("Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of", edge_dist, "between embedded values\n"),
                       no = paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n"))
size_nodes_by = count_col # can use a double, e.g., 1.0, for fixed size
node_size_limits = NULL # numeric, length 2
color_nodes_by = NULL # use NULL to automatically determine
color_scheme = "default"

# Output Settings
output_dir = file.path(dir_main, "output") # if NULL, output is not saved to file
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




#### INPUT CHECKS ####

# Atchley factor embedding only applicable to amino acid sequences
if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
  stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }

# aggregate_reads only applicable to amino acid sequences


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
  cluster_id_col <- degree_col <- NULL
  if ("cluster_id" %in% names(data)) { cluster_id_col <- "cluster_id" }
  if ("degree" %in% names(data)) { degree_col <- "degree" }
  cluster_info <- getClusterStats(
    data, adjacency_matrix, clone_seq_col, count_col,
    cluster_id_col, degree_col) }

clone_col = clone_seq_col # name or number of column of `data` containing the clone sequences

cluster_id_col <- degree_col <- NULL
if ("cluster_id" %in% names(data)) { cluster_id_col <- "cluster_id" }
if ("degree" %in% names(data)) { degree_col <- "degree" }


seq_length_col = NULL




# Function ----------------------------------------------------------------

cluster_info <- getClusterStats(
  data, adjacency_matrix, clone_seq_col, count_col,
  cluster_id_col,
  degree_col)
