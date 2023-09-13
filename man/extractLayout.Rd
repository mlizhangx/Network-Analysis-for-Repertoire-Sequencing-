\name{extractLayout}
\alias{extractLayout}
\title{
  Get Coordinate Layout From Graph Plot
}
\description{
  Given a \code{\link[ggraph]{ggraph}} plot, extract the coordinate layout of
  the graph nodes as a two-column matrix.
}
\usage{
extractLayout(plot)
}
\arguments{
  \item{plot}{An object of class \code{\link[ggraph]{ggraph}}.
}
}
\details{
  Equivalent to \code{as.matrix(plot$data[c("x", "y")])}.
}
\value{
  A matrix with two columns and one row per network node. Each row contains the
  Cartesian coordinates of the corresponding node.
}
\references{
Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang. (2023).
NAIR: Network Analysis of Immune Repertoire. \emph{Frontiers in Immunology}, vol. 14.
\href{https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full}{doi: 10.3389/fimmu.2023.1181825}

\href{https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/index.html}{Webpage for the NAIR package}
}
\author{
    Brian Neal (\email{Brian.Neal@ucsf.edu})
}
\examples{

set.seed(42)
toy_data <- simulateToyData()
net <- buildRepSeqNetwork(toy_data, "CloneSeq", print_plots = TRUE)

my_layout <- extractLayout(net$plots[[1]])

# same as `graph_layout` element in the plot list
all.equal(my_layout, net$plots$graph_layout, check.attributes = FALSE)

}
