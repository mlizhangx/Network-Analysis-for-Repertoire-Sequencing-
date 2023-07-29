\name{saveNetworkPlots}
\alias{saveNetworkPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Write Plots to a PDF
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given a list of plots, write all plots to a single pdf file.
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
    A list whose elements can be written to pdf using the \code{\link[grDevices]{pdf}} function, such as the \code{\link[ggplot2]{ggplot}} class.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
    The file path for the pdf.
}
  \item{pdf_width}{
%%     ~~Describe \code{plot_width} here~~
    Passed to the \code{width} argument of the \code{\link[grDevices]{pdf}} function.
}
  \item{pdf_height}{
%%     ~~Describe \code{plot_height} here~~
    Passed to the \code{height} argument of the \code{\link[grDevices]{pdf}} function.
}
}
\author{
%%  ~~who you are~~
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{
\dontrun{
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
