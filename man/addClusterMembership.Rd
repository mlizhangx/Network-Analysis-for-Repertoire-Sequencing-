\name{addClusterMembership}
\alias{addClusterMembership}
\title{
    Partition a Network Graph Into Clusters
}
\description{
    Given a list of network objects returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
    \code{\link[=generateNetworkObjects]{generateNetworkObjects()}},
    partitions the network graph into clusters using the specified clustering
    algorithm, adding a cluster membership variable to the node metadata.
}
\usage{
addClusterMembership(
  net,
  cluster_fun = "fast_greedy",
  cluster_id_name = "cluster_id",
  overwrite = FALSE,
  verbose = FALSE,
  ...,
  data = deprecated(),
  fun = deprecated()
)
}
\arguments{
  \item{net}{
      A \code{\link{list}} of network objects conforming to the output of
      \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} or
      \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}.
      See details.
      Alternatively, this argument accepts the network
      \code{\link[igraph]{igraph}}, with the node metadata passed to the
      \code{data} argument. However, this alternative functionality is
      deprecated and will eventually be removed.
}
  \item{cluster_fun}{
      A character string specifying the clustering algorithm to use.
      See details.
}
  \item{cluster_id_name}{
      A character string specifying the name of the cluster membership variable
      to be added to the node metadata.
}
  \item{overwrite}{
      Logical. Should the variable specified by \code{cluster_id_name} be
      overwritten if it already exists?
  }
  \item{verbose}{
      Logical. If \code{TRUE}, generates messages about the tasks
      performed and their progress, as well as relevant properties of intermediate
      outputs. Messages are sent to \code{\link[=stderr]{stderr()}}.
  }
  \item{...}{
      Named optional arguments to the function specified by \code{cluster_fun}.
  }
  \item{data}{
      \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
      See \code{net}.
}
  \item{fun}{
      \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
      Replaced by \code{cluster_fun}.
}
}
\details{
    The list \code{net} must contain the named elements
    \code{igraph} (of class \code{\link[igraph]{igraph}}),
    \code{adjacency_matrix} (a \code{\link{matrix}} or
    \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}} encoding edge connections),
    and \code{node_data} (a \code{\link{data.frame}} containing node metadata),
    all corresponding to the same network. The lists returned by
    \code{\link[=buildRepSeqNetwork]{buildRepSeqNetwork()}} and
    \code{\link[=generateNetworkObjects]{generateNetworkObjects()}}
    are examples of valid inputs for the \code{net} argument.

    Alternatively, the \code{igraph} may be passed to \code{net} and the node
    metadata to \code{data}. However, this alternative functionality is
    deprecated and will eventually be removed.

    A clustering algorithm is used to partition the network graph into clusters
    (densely-connected subgraphs). Each cluster represents a collection of
    clones/cells with similar receptor sequences. The method used to partition
    the graph depends on the choice of clustering algorithm, which is specified
    using the \code{cluster_fun} argument.

    The available options for \code{cluster_fun} are listed below. Each refers
    to an \code{\link[igraph:igraph-package]{igraph}} function implementing a
    particular clustering algorithm. Follow the links to learn more about the
    individual clustering algorithms.

    \itemize{
    \item \code{\link[igraph:cluster_edge_betweenness]{"edge_betweenness"}}
    \item \code{\link[igraph:cluster_fast_greedy]{"fast_greedy"}}
    \item \code{\link[igraph:cluster_infomap]{"infomap"}}
    \item \code{\link[igraph:cluster_label_prop]{"label_prop"}}
    \item \code{\link[igraph:cluster_leading_eigen]{"leading_eigen"}}
    \item \code{\link[igraph:cluster_leiden]{"leiden"}}
    \item \code{\link[igraph:cluster_louvain]{"louvain"}}
    \item \code{\link[igraph:cluster_optimal]{"optimal"}}
    \item \code{\link[igraph:cluster_spinglass]{"spinglass"}}
    \item \code{\link[igraph:cluster_walktrap]{"walktrap"}}
    }

    Optional arguments to each clustering algorithm can have their
    values specified using the ellipses (\code{...}) argument of
    \code{addClusterMembership()}.

    Each cluster is assigned a numeric cluster ID. A cluster membership variable,
    whose name is specified by \code{cluster_id_name}, is added to the node
    metadata, encoding the cluster membership of the node for each row. The cluster
    membership is encoded as the cluster ID number of the cluster to which the node
    belongs.

    The \code{overwrite} argument controls whether to overwrite pre-existing data.
    If the variable specified by \code{cluster_id_name} is already present in
    the node metadata, then \code{overwrite} must be set to \code{TRUE} in
    order to perform clustering and overwrite the variable with new cluster
    membership values. Alternatively, by specifying a value for
    \code{cluster_id_name} that is not among the variables in the node metadata,
    a new cluster membership variable can be created while preserving the old
    cluster membership variable. In this manner, clustering can be performed
    multiple times on the same network using different clustering algorithms,
    without losing the results.
}

\value{
    If the variable specified by \code{cluster_id_name} is not present in
    \code{net$node_data}, returns a copy of \code{net} with this variable
    added to \code{net$node_data} encoding the cluster membership of the network
    node corresponding to each row. If the variable is already present and
    \code{overwrite = TRUE}, then its values are replaced with the new values
    for cluster membership.

    Additionally, if \code{net} contains a list named \code{details}, then the
    following elements will be added to \code{net$details} if they do not
    already exist:

    \item{\code{clusters_in_network}}{
        A named numeric vector of length 1. The first entry's name is the
        name of the clustering algorithm, and its value is the number of clusters
        resulting from performing clustering on the network.
    }
    \item{\code{cluster_id_variable}}{
        A named numeric vector of length 1. The first entry's name is the
        name of the clustering algorithm, and its value is the name of the
        corresponding cluster membership variable in the node metadata
        (i.e., the value of \code{cluster_id_name}).
    }

    If \code{net$details} already contains these elements, they will be updated
    according to whether the cluster membership variable specified by
    \code{cluster_id_name}
    is added to \code{net$node_data} or already exists and is overwritten.
    In the former case (the cluster membership variable does not already exist),
    the length of each vector
    (\code{clusters_in_network}) and (\code{cluster_id_variable})
    is increased by 1, with the new information appended as a new named entry
    to each. In the latter case (the cluster membership variable is overwritten),
    the new information overwrites the name and value of the last entry of each
    vector.

    In the event where \code{overwrite = FALSE} and \code{net$node_data} contains
    a variable with the same name as the value of \code{cluster_id_name}, then an
    unaltered copy of \code{net} is returned with a message notifying the user.

    Under the alternative (deprecated) input format where the node metadata is
    passed to \code{data} and the \code{igraph} is passed to \code{net}, the
    node metadata is returned instead of the list of network objects, with the
    cluster membership variable added or updated as described above.
}
\seealso{
\code{\link[=addClusterStats]{addClusterStats()}}
\code{\link[=addClusterLabels]{addClusterLabels()}}
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

net <- generateNetworkObjects(
  toy_data, "CloneSeq"
)

# Perform cluster analysis,
# add cluster membership to net$node_data
net <- addClusterMembership(net)

net$details$clusters_in_network
net$details$cluster_id_variable

# overwrite values in net$node_data$cluster_id
# with cluster membership values obtained using "cluster_leiden" algorithm
net <- addClusterMembership(
  net,
  cluster_fun = "leiden",
  overwrite = TRUE
)

net$details$clusters_in_network
net$details$cluster_id_variable

# perform clustering using "cluster_louvain" algorithm
# saves cluster membership values to net$node_data$cluster_id_louvain
# (net$node_data$cluster_id retains membership values from "cluster_leiden")
net <- addClusterMembership(
  net,
  cluster_fun = "louvain",
  cluster_id_name = "cluster_id_louvain",
)

net$details$clusters_in_network
net$details$cluster_id_variable

}
