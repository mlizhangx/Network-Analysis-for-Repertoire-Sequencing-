\name{addGraphLabels}
\alias{addGraphLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
  Label Nodes on a Graph Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Annotate a graph plot to add custom labels to the nodes.
}
\usage{
addGraphLabels(plot, node_labels, size = 5, color = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plot}{
%%     ~~Describe \code{plot} here~~
  A ggraph object.
}
  \item{node_labels}{
%%     ~~Describe \code{node_labels} here~~
  A vector of node labels (the length should match the number of nodes in the plot).
}
  \item{size}{
%%     ~~Describe \code{size} here~~
  The size of the node labels.
}
  \item{color}{
%%     ~~Describe \code{color} here~~
  The color of the node labels.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
  Labels are added via \code{ggraph::geom_node_text}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  A \code{ggraph} object; the input plot with labels added.
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
%% ~Make other sections like Warning with \section{Warning }{....} ~
