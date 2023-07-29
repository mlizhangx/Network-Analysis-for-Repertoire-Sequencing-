\name{saveNetwork}
\alias{saveNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
    Save List of Network Objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    Given a list of network objects such as that returned by \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}, saves
    its contents according to the specified file format scheme.
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
    A list conforming to the format returned by \code{\link{buildRepSeqNetwork}} or \code{\link{generateNetworkObjects}}. See details.
}
  \item{output_dir}{
%%     ~~Describe \code{output_dir} here~~
    A file path specifying the directory in which to save the file(s).
}
  \item{output_type}{
%%     ~~Describe \code{output_dir} here~~
  A character string specifying the file format to use when writing output to file. Valid options are \code{"individual"}, \code{"rds"} and \code{"rda"}.
}
  \item{output_filename}{
%%     ~~Describe \code{output_filename} here~~
  A character string to be used as a common filename prefix for any files saved.
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
\details{
%%  ~~ If necessary, more details than the description above ~~
    The list provided to the \code{net} argument must contain at least the components \code{igraph} (an igraph object), \code{adjacency_matrix} (a matrix or dgCMatrix) and \code{node_data} (a data frame); the components \code{cluster_data} (a data frame) and \code{plots} (a list of \code{ggraph} objects) are optional, and will also be saved if present.

    The  saves each item as a separate, uncompressed file, with data frames saved in csv format. Setting \code{output_type} to \code{"rda"} or \code{"rds"} saves the output list as a rda or rds file, respectively. For all output types, any plots present will also be saved in a pdf file.

    Under the default file format scheme \code{output_type = "individual"}, the igraph is saved in edgelist
    format as a .txt file; the adjacency matrix is saved either in .csv format
    (for dense matrices) or .mtx format (for sparse matrices); and plots are
    saved only in a .pdf file, meaning that any plots will need to be generated
    again in order to be modified. If \code{output_type} is \code{"rds"} or
    \code{"rda"}, the entire list \code{net} will be saved (the R symbol for the
    object will be \code{net} if saving as .rda).
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

# Save network objects
saveNetwork(net, output_dir = tempdir())
}
}