\name{saveNetworkPlots}
\alias{saveNetworkPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Save List of GGPlots to a Single PDF
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Print a list of \code{ggplot} objects to a single pdf file.
}
\usage{
saveNetworkPlots(plotlist,
                 outfile = "MyRepSeqNetwork.pdf",
                 pdf_width = 12, pdf_height = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plotlist}{
%%     ~~Describe \code{plotlist} here~~
    A list of \code{ggplot} objects.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
    The file to which the printed plots should be written.
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
%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
%% ~put references to the literature/web site here ~
%%    Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang. "Network Analysis of Immune Repertoire (NAIR) with Advanced Machine Learning Techniques." In: Briefings in Bioinformatics (Submitted: July 18, 2022).

    https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-
}
\author{
%%  ~~who you are~~
    Brian Neal, Hai Yang, Jason Cham, Zenghua Fan, Tao He and Li Zhang.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{plotNetworkGraph}}
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

# Generate network for data
net <- generateNetworkObjects(data, "clone_seq")

# Plot network graph
net$plots <- generateNetworkGraphPlots(
  net$igraph, net$node_data, color_nodes_by = c("group", "count"))

# Save plots
saveNetworkPlots(net$plots, outfile = file.path(tempdir(), "network.pdf"))
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
