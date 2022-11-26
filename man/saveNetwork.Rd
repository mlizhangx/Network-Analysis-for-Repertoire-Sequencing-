\name{saveNetwork}
\alias{saveNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Save List of RepSeq Network Objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given a list conforming to the output of \code{buildRepSeqNetwork()}, saves
    its contents according to the specified file format scheme. Also saves
    any graph plots present to a single pdf.
}
\usage{
saveNetwork(net, output_dir = getwd(),
            output_type = "individual",
            output_filename = "MyRepSeqNetwork",
            pdf_width = 12, pdf_height = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
%%     ~~Describe \code{net} here~~
    A list conforming to the format of the value returned by \code{buildRepSeqNetwork}. Must contain at least the components \code{igraph} (an igraph object), \code{adjacency_matrix} (a matrix or dgCMatrix) and \code{node_data} (a data frame); the components \code{cluster_data} (a data frame) and \code{plots} (a list of \code{ggraph} objects) are optional, and will also be saved if present.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    The directory in which to save the file(s).
}
  \item{output_type}{
%%     ~~Describe \code{output_dir} here~~
  A character string specifying the file format to use when writing output to file. Default \code{"individual"} saves each item as a separate, uncompressed file, with data frames saved in csv format. \code{"rda"} and \code{"rds"} save the output list as a rda and rds file, respectively. For all output types, any plots present will also be saved in a pdf file.
}
  \item{output_filename}{
%%     ~~Describe \code{output_filename} here~~
  A character string to be used as a common filename prefix for any files saved.
}
  \item{pdf_width}{
%%     ~~Describe \code{plot_width} here~~
    Passed to the \code{width} argument of \code{grDevices::pdf()}.
}
  \item{pdf_height}{
%%     ~~Describe \code{plot_height} here~~
    Passed to the \code{height} argument of \code{grDevices::pdf()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
    When \code{output_type = "individual"}, the igraph is saved in edgelist
    format as a .txt file; the adjacency matrix is saved either in .csv format
    (for dense matrices) or .mtx format (for sparse matrices); and plots are
    saved only in a .pdf file, meaning that any plots will need to be generated
    again in order to be modified. If \code{output_type} is \code{"rds"} or
    \code{"rda"}, the entire list \code{net} will be saved (the R symbol for the
    object will be \code{net} if saving as .rda).
}
\references{
%% ~put references to the literature/web site here ~
%%    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

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
data <- data.frame(seqs = clone_seqs)
net <- generateNetworkObjects(data, "seqs")
saveNetwork(net, output_dir = tempdir())

}