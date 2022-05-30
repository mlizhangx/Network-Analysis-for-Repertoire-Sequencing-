
# Setup working environment -----------------------------------------------

library(RepSeqNetworkAnalysis)

# set working directory to dir of current script
dir_main <- dirname(rstudioapi::getActiveDocumentContext()$path)

# load data sets
data <- read.table(
  file.path(dir_main, "TRB-Pt-10-2-500ng-15-04-2020-gDNA_S48.clones.txt"),
  sep = '\t', header = TRUE)

setwd(file.path(dir_main, "output"))


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
cluster_stats = FALSE
node_stat_settings = node_stat_settings(cluster_id = cluster_stats)

# Plot Settings
plot_title = paste("RepSeq network on CDR3",
                   clone_seq_type, "sequence")
plot_subtitle = ifelse(dist_type == "euclidean_on_atchley",
                       yes = paste("Clone sequences embedded in Euclidean 30-space based on Atchley factor representation using deep learning\nEdges based on a maximum Euclidean distance of", edge_dist, "between embedded values\n"),
                       no = paste("Edges based on a maximum", dist_type, "distance of", edge_dist, "\n"))
size_nodes_by = count_col # can use a double, e.g., 1.0, for fixed size
node_size_limits = NULL # numeric, length 2
color_nodes_by = NULL # use NULL to automatically determine
color_scheme = "default"

# Output Settings
output_dir = getwd() # if NULL, output is not saved to file
plot_outfile = "network_graph.pdf"
data_outfile = "node_data.csv"
igraph_outfile = NULL # .txt
matrix_outfile = NULL # .mtx (.csv for euclidean on atchley)
return_all = FALSE # if false, only the node data is returned (unless cluster_stats = TRUE and output_dir = NULL, in which case a list containing the node data and cluster info is returned, with a warning)


# Atchley factor embedding only applicable to amino acid sequences
if (dist_type == "euclidean_on_atchley" & clone_seq_type != "amino acid") {
  stop("distance type 'euclidean_on_atchley' only applicable to amino acid sequences") }


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




# Function internals ------------------------------------------------------


# Convert column specifications from numeric to character if not already
if (is.numeric(clone_col)) { clone_col <- names(data)[clone_col] }
if (is.numeric(count_col)) { count_col <- names(data)[count_col] }
if (is.numeric(freq_col)) { freq_col <- names(data)[freq_col] }

# Define grouping variable(s)
grouping_variables <- list(data[ , clone_col])
names(grouping_variables) <- clone_col
if (!is.null(grouping_cols)) {
  if (is.numeric(grouping_cols)) { grouping_cols <- names(data)[grouping_cols] }
  for (i in 1:length(grouping_cols)) {
    grouping_variables$newvar <- data[ , grouping_cols[[i]]]
    names(grouping_variables)[[length(grouping_variables)]] <-
      grouping_cols[[i]] } }

# aggregate the reads by group
data_to_aggregate <- list("aggCloneCount" = data[ , c(count_col)],
                          "aggCloneFreq" = data[ , c(freq_col)])
agg_counts <- stats::aggregate(data_to_aggregate,
                               by = grouping_variables, FUN = sum)

# add variable for num reads (row count)
groups <- as.data.frame(grouping_variables)
names(groups)[[1]] <- "temporary_placeholder_name" # for summarize function
num_reads <- dplyr::summarize(dplyr::group_by_all(groups),
                              numReads = length(temporary_placeholder_name))
names(num_reads)[[1]] <- clone_col # replace placeholder name with orig

# Merge aggregate counts with num reads
out <- merge(agg_counts, num_reads, by = c(clone_col, grouping_cols))
