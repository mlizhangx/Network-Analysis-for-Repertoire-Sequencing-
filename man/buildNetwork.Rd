\name{buildNetwork}
\alias{buildNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate Network Graph and Statistics From RepSeq Data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given clonotype sequences and counts from repertoire-sequence data, builds the corresponding network graph using `igraph`, computes node-level and cluster-level network statistics, and creates a `ggraph` visual plot of the network.
}
\usage{
%- maybe also 'usage' for other objects documented here.
buildNetwork(clonotypes, counts, frequencies,
             group_labels = NULL,
             distance_type = "levenshtein",
             clone_seq_type = "unspecified clonotype",
             sample_name = "unknown_sample",
             group_type = "no group labels provided",
             save_output = FALSE,
             report_time = FALSE)
}
\arguments{
  \item{clonotypes}{
%%     ~~Describe \code{clonotypes} here~~
A character vector of clonotype sequences (e.g., nucleotide or amino acid sequences), where duplicate sequences are allowed.
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
A numeric vector of the same length as `clonotypes`, containing the corresponding clonotype counts (e.g., of templates/reads).
}
  \item{frequencies}{
%%     ~~Describe \code{frequencies} here~~
A numeric vector of the same length as `clonotypes`, containing the corresponding relative frequencies.
}
  \item{group_labels}{
%%     ~~Describe \code{group_labels} here~~
An optional vector of the same length as `clonotypes`, representing a categorical variable of interest (e.g., treatment vs. control status, before vs. after treatment, etc.). If supplied, a variable 'groupLabel' will be included in the node-level data, reporting the level of `group_labels` with the greatest aggregate count for the clonotype sequence corresponding to each node; nodes in the `ggraph` will be colored according to this variable.
}
  \item{distance_type}{
%%     ~~Describe \code{distance_type} here~~
An optional character string specifying the dissimilarity measure to use when building the network. Default is "levenshtein" for Levenshtein distance. The other option currently supported is the Hamming distance, which can be used with `distance_type = "hamming"`. The Hamming distance assumes clonotype sequences of equal length; sequences shorter than the maximum length will have space characters appended to pad the missing length.
}
  \item{clone_seq_type}{
%%     ~~Describe \code{clone_seq_type} here~~
An optional character string specifying the type of clonotype sequence used (e.g., nucleotide, amino acid).  This has no effect other than altering the description of the output and the label used in the `ggraph` plot.
}
  \item{sample_name}{
%%     ~~Describe \code{sample_name} here~~
An optional character string specifying the name or designation of the sample of repertoire-sequence data.  This has no effect other than altering the description of the output and the label used in the `ggraph` plot.
}
  \item{group_type}{
%%     ~~Describe \code{group_type} here~~
An optional character string specifying the name of the categorical variable supplied in the `group_labels` argument.  This has no effect other than altering the description of the output and the label used in the `ggraph` plot, and only applies if `group_labels` is non-null.
}
  \item{save_output}{
%%     ~~Describe \code{save_output} here~~
An optional logical expression evaluating to 'TRUE' or 'FALSE', indicating whether the output should be automatically saved. If 'TRUE', the network graph is saved in edgelist format as a text file, the node-level and cluster-level data frames are saved as csv files, and the `ggraph` is plotted and saved as a pdf. A subdirectory of the current working directory is automatically created to contain the saved files, distinguished according to `sample_name`, `clone_seq_type` and `distance_type`.
}
  \item{report_time}{
%%     ~~Describe \code{report_time} here~~
An optional logical expression evaluating to 'TRUE' or 'FALSE', indicating whether to report the time taken to process the data, build the network graph and compute the network statistics.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Before building the network, the function first cleans and aggregates the clonotype sequences and counts from the repertoire-sequence data, removing sequences with fewer than 3 characters. The function then computes the pairwise dissimilarity between the clonotype sequences and constructs a corresponding adjacency matrix for the sequences corresponding to nodes with positive network degree. An igraph is generated using the adjacency matrix, and the node-level and cluster-level network statistics are computed and stored in separate data frames. A ggraph is created, which can be printed to display a visual plot of the network graph; optionally, the nodes can be colored according to a categorical variable supplied via the `group_labels` argument (the level of `group_labels` assigned to each node is the value with the greatest aggregate count for the corresponding clonotype sequence).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing the following items:
\item{settings}{A list of parameter values used when calling the `buildNetwork()` function.}
\item{igraph}{An `igraph` network object containing the edge list for the network.}
\item{cell_info}{A data frame containing the node-level network characteristics, along with the counts and frequencies.}
\item{cluster_info}{A data frame containing the cluster-level network characteristics.}
\item{ggraph}{A `ggraph` object containing the visual plot of the network graph.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Li Zhang, Hai Yang, Brian Neal
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{


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
