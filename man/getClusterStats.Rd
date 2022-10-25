\name{getClusterStats}
\alias{getClusterStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Compute Cluster-Level Network Statistics
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Given immune repertoire sequence data and the adjacency matrix for the
  corresponding repertoire network, perform clustering and compute cluster-level
  network statistics.
}
\usage{
getClusterStats(data,
                adjacency_matrix,
                seq_col,
                count_col = NULL,
                cluster_id_col = NULL,
                degree_col = NULL,
                seq_length_col = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
      A data frame containing the immune repertoire sequence data, with clones
      indexed by row.
}
  \item{adjacency_matrix}{
%%     ~~Describe \code{adjacency_matrix} here~~
    A square adjacency matrix corresponding to the repertoire network for the clones
    in `data`, with row dimension matching that of `data`. Sparse matrices
    created using the `Matrix` package are supported.
}
  \item{seq_col}{
%%     ~~Describe \code{seq_col} here~~
      The column name or number of `data` containing the receptor sequences (e.g.,
      TCR CDR3 nucleotide sequence or amino acid sequence).
}
  \item{count_col}{
%%     ~~Describe \code{count_col} here~~
      Optional. The column name or number of `data` containing the clone counts. If supplied, then related cluster-level properties will be computed.
}
  \item{cluster_id_col}{
%%     ~~Describe \code{cluster_id_col} here~~
      The column name or number of `data` containing the cluster membership for each clone. If `NULL`, clustering will be performed using `igraph::cluster_fast_greedy()` in order to determine cluster membership.
}
  \item{degree_col}{
%%     ~~Describe \code{degree_col} here~~
      The column name or number of `data` containing the network degree for each clone. If `NULL`, the network degree will be computed using `igraph::degree`.
}
  \item{seq_length_col}{
%%     ~~Describe \code{seq_length_col} here~~
      The column name or number of `data` containing the sequence length for each clone. If `NULL`, the values will be computed based on the sequences in `seq_col`.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    The cluster-level statistics `diameter_length`, `assortativity`, `transitivity`,
    `edge_density`, `degree_centrality_index`, `closeness_centrality_index`,
    `eigen_centrality_index`, and `eigen_centrality_eigenvalue` are computed for each cluster
    by first generating the network graph for the cluster using the adjacency
    submatrix corresponding to that cluster. As a result, computing cluster statistics
    can potentially take some time for very large networks with many large clusters.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    A data frame containing one row for each cluster and the following columns:

  \item{cluster_id}{The cluster ID number.}
  \item{node_count}{The number of nodes in the cluster.}
  \item{mean_seq_length}{The mean sequence length in the cluster.}
  \item{mean_degree}{The mean network degree in the cluster.}
  \item{max_degree}{The maximum network degree in the cluster.}
  \item{seq_w_max_degree}{The clone sequence possessing the maximum degree within the cluster.}
  \item{agg_clone_count}{The aggregate clone count for all clones in the cluster.}
  \item{max_clone_count}{The maximum clone count among all clones in the cluster.}
  \item{seq_w_max_count}{The clone sequence possessing the maximum clone count within the cluster.}
  \item{diameter_length}{The longest geodesic distance in the cluster, computed using `length(igraph::get_diameter())`.}
  \item{assortativity}{The assortativity coefficient of the cluster's graph, based on the degree (minus one) of each node in the cluster (with the degree computed based only upon the nodes within the cluster). Computed using `igraph::assortativity_degree`.}
  \item{transitivity}{The transitivity (i.e., clustering coefficient) for the cluster's graph, which estimates the probability that adjacent vertices are connected. Computed using `igraph::transitivty()` with `type = "global"`.}
  \item{edge_density}{The number of edges in the cluster as a fraction of the maximum possible number of edges. Computed using `igraph::edge_density()`.}
  \item{degree_centrality_index}{The cluster-level centrality index based on degree within the cluster graph. Computed using `igraph::centr_degree()$centralization`.}
  \item{closeness_centrality_index}{The cluster-level centrality index based on closeness, i.e., distance to other nodes in the cluster. Computed using `igraph::centr_clo()$centralization`.}
  \item{eigen_centrality_index}{The cluster-level centrality index based on the eigenvector centrality scores, i.e., values of the first eigenvector of the adjacency matrix for the cluster. Computed using `igraph::centr_eigen()$centralization`.}
  \item{eigen_centrality_eigenvalue}{The eigenvalue corresponding to the first eigenvector of the adjacency matrix for the cluster. Computed using `igraph::eigen_centrality()$value`.}
}
\references{
%% ~put references to the literature/web site here ~
    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\examples{

# Generate some data
set.seed(42)
sample_size <- 200
group_labels <- rep("Control", times = sample_size)
assign_group <- sample(c(TRUE, FALSE), size = sample_size, replace = TRUE)
base_seq_length <- 7
characters <- c("G", "A", "T", "C")
char_probs <- c(2/3, 2/9, 1/18, 1/18)
characters_sample <- sample(characters,
                            size = base_seq_length * sample_size,
                            replace = TRUE, prob = char_probs)
characters_sample <- matrix(characters_sample, ncol = sample_size)
clone_seqs <- apply(characters_sample, MARGIN = 2,
                    FUN = function(x) paste0(x, collapse = ""))
clone_seqs_append <- rep("AATC", times = sample_size)
append_latent_prob <- runif(sample_size, min = 0, max = 1)
for (i in 1:sample_size) {
  case_group <- assign_group[[i]]
  if (case_group) group_labels[[i]] <- "Case"
  if ((case_group & append_latent_prob[[i]] > 0.9) |
      (!case_group & append_latent_prob[[i]] > 0.5)) {
    clone_seqs_append[[i]] <- "AATCGGGG"
  } else if ((case_group & append_latent_prob[[i]] > 0.8) |
             (!case_group & append_latent_prob[[i]] > 0.3)) {
    clone_seqs_append[[i]] <- "AATCGGT"
  } else if ((case_group & append_latent_prob[[i]] > 0.7) |
             (!case_group & append_latent_prob[[i]] > 0.2)) {
    clone_seqs_append[[i]] <- "AATCGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.5) |
             (!case_group & append_latent_prob[[i]] > 0.1)) {
    clone_seqs_append[[i]] <- "AATTGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.25) |
             (!case_group & append_latent_prob[[i]] > 0.05)) {
    clone_seqs_append[[i]] <- "AATTG"
  }
  clone_seqs[[i]] <- paste0(clone_seqs[[i]], clone_seqs_append[[i]],
                            collapse = "")
}
counts <- rbinom(sample_size, size = 300, prob = 0.1)
frequencies <- counts/sum(counts)
data <- data.frame("clone_seq" = clone_seqs,
                   "count" = counts,
                   "frequency" = frequencies,
                   "group" = group_labels)

# Generate adjacency matrix for data
adjmat <-
  generateNetworkFromSeqs(
    data$clone_seq,
    drop_isolated_nodes = FALSE,
    return_type = "adjacency_matrix"
  )

# Get cluster stats
cluster_info <-
  getClusterStats(data, adjmat,
                  seq_col = "clone_seq",
                  count_col = "count")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
