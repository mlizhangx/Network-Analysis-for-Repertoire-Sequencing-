\name{generateNetworkFromSeqs}
\alias{generateNetworkFromSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Generate the Network for a List of Receptor Sequences
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given a list of receptor sequences for an immune repertoire, computes the
    adjacency matrix for the repetroire network based on sequence similarity,
    then generates the corresponding network.
}
\usage{
generateNetworkFromSeqs(
    seqs,
    dist_type = "hamming",
    dist_cutoff = 1,
    drop_isolated_nodes = TRUE,
    contig_ids = seq_along(seqs),
    outfile_adjacency_matrix = NULL,
    outfile_distance_matrix = NULL,
    return_type = "network")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
%%     ~~Describe \code{seqs} here~~
    A character vector containing the list of receptor sequences.
}
  \item{dist_type}{
%%     ~~Describe \code{dist_type} here~~
    The distance metric to use as a measure of similarity between two receptor
    sequences. Valid options are \code{"hamming"} (the default), \code{"levenshtein"}
    (i.e., edit distance) and \code{"euclidean_on_atchley"} (sequences are encoded
    as numeric vectors using deep learning, and the Euclidean distance is measured
    between the encoded vectors; only applicable to TCR CDR3 amino acid sequences).
}
  \item{dist_cutoff}{
%%     ~~Describe \code{edge_dist} here~~
      Numeric; the maximum distance allowed between two receptor sequences that share a network edge.
}
  \item{drop_isolated_nodes}{
%%     ~~Describe \code{drop_isolated_nodes} here~~
      Logical. If TRUE, nodes that are not incidental to any network edges (i.e., cells/clones whose receptor sequences are not similar to any others) are dropped from the network.
}
  \item{contig_ids}{
%%     ~~Describe \code{contig_ids} here~~
      Passed to `embedTCRSeqsByAtchleyFactor()`.
}
  \item{outfile_adjacency_matrix}{
%%     ~~Describe \code{outfile_adjacency_matrix} here~~
      The name of the file where the adjacency matrix will be saved. This should end with the `.mtx` file extension. If `NULL` (default), the adjacency matrix will not be saved. Only applicable when `dist_type` is `"hamming"` or `"levenshtein"`.
}
  \item{outfile_distance_matrix}{
%%     ~~Describe \code{outfile_distance_matrix} here~~
      The name of the file where the matrix of pairwise distances between sequences will be saved. This should end with the `.csv` file extension. If `NULL` (default), the distance matrix will not be saved. Only applicable when `dist_type = "euclidean_on_atchley"`, as explicit distances are not necessarily computed for other distance types.
}
  \item{return_type}{
%%     ~~Describe \code{return_type} here~~
    Character string specifying the type of output. Valid options are `"network"` (the default), which returns an `igraph` network object, and `"adjacency_matrix"`, which returns the adjacency matrix.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
    If `return_type = "network"`, an `igraph` object; if `return_type = "adjacency_matrix"`,
    then a matrix whose class depends on `dist_type`: If `dist_type = "euclidean_on_atchley"`,
    a dense matrix of type `"double"`, otherwise a sparse matrix with Formal class `"dgCMatrix"`.
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

net <- generateNetworkFromSeqs(seqs = clone_seqs,
                               dist_type = "levenshtein")
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
