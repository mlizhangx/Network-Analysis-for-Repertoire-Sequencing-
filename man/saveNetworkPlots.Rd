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
# Simulate some toy data for demonstration
toy_data <- simulateToyData()

# Generate network for data
net <- generateNetworkObjects(toy_data, "CloneSeq")

# Plot network graph
net$plots <- generateNetworkGraphPlots(
  net$igraph, net$node_data, color_nodes_by = c("SampleID", "CloneCount"))

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
