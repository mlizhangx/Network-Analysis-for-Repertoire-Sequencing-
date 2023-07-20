
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

We introduce the `NAIR` package, a powerful tool for analyzing the
adaptive immune repertoire using network analysis based on similarities
among receptor sequences, which implements methods from the following
paper:

Hai Yang, Jason Cham, Zenghua Fan, Brian Neal, Tao He and Li Zhang.
NAIR: Network Analysis of Immune Repertoire. *Frontiers in Immunology*
(forthcoming).

Our immune repertoire analysis tool allows users to perform various
network analysis tasks on AIRR-Seq data. This includes computing local
and global network properties of nodes and clusters, which can provide
insights into the structural organization of the immune repertoire
network.

In addition, our tool enables users to search across multiple AIRR-Seq
samples for clones/clusters associated with subject characteristics,
disease conditions or clinical outcomes, as well as identify public
clones/clusters. This can help researchers identify potentially
important TCR/BCR clones.

To aid in visualization and interpretation of the immune repertoire
network, our tool also allows for customized visualizations of the
network. Furthermore, our tool enables users to perform further
downstream analysis on the immune repertoire data, such as clustering
analysis and differential abundance testing. These downstream analyses
can provide additional insights into the immune response.

### What data does `NAIR` support?

`NAIR` supports bulk and single-cell immune repertoire sequence data for
T-cell or B-cell receptors (TCR or BCR).

- **Single-cell data:** Each row is a single cell
- **Bulk data:** Each row is a distinct TCR/BCR clone (unique
  combination of V-D-J genes and nucleotide sequence) and typically
  includes a corresponding measurement of clonal abundance (e.g., clone
  count and clone frequency/fraction)

### How does `NAIR` model the immune repertoire as a network?

- Each cell (single-cell data) or clone (bulk data) is modeled as a node
  (vertex) in the network
- For each node, we consider the corresponding receptor sequence
  (nucleotide or amino acid)
- For each pair of nodes, we measure the similarity in their receptor
  sequences (using the Hamming or Levenshtein distance)
- An edge is drawn between two nodes if the distance is below a
  specified threshold
  - For single-cell data, sequences from two chains (e.g., alpha chain
    and beta chain) can be jointly used to determine similarity between
    cells, considering cells as similar when the sequences for both
    chains are similar (i.e., when the distance for each chain is below
    the threshold)
- Clustering analysis is used to partition the network graph into
  clusters (densely-connected subgraphs)
  - Many clustering algorithms are available, with each seeking to
    identify the “best” configuration of clusters according to different
    graph criteria
- Network statistics characterize the repertoire in terms of the local
  and global structural properties of its graph
- Customized visual plots of the network graph are generated, with nodes
  colored according to desired meta-data (e.g., disease status, sample
  ID, cluster ID, clone count, etc.)

# Installation

The current development version of `NAIR` can be installed from github
using the following commands:

``` r
install.packages("devtools")
devtools::install_github(
  "mlizhangx/Network-Analysis-for-Repertoire-Sequencing-",
  build_vignettes = TRUE)
```

Installing the development version requires a toolchain compiler. On
Windows, this means downloading and installing Rtools. On MacOS, this
entails installing XCode Command Line Tools (“XCode CLI”) and the
correct version of gfortran for your macOS version (instructions
[here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/)).

# Documentation

Once the package is installed, type the following commands to access the
package vignettes and documentation:

``` r
# Package introduction vignette
vignette("NAIR", package = "NAIR")

# Display vignettes in html browser
browseVignettes("NAIR")

# Browse documentation
help(package = "NAIR")
```

# The `buildRepSeqNetwork` function

General network analysis on AIRR-Seq data is performed using the
`buildRepSeqNetwork` function. This function does the following:

- Filters the input AIRR-Seq data according to user specifications
- Builds the network graph for the immune repertoire
- Performs desired network analysis
  - Cluster analysis
  - Computation of network properties
- Generates customized visual plot of the network graph using `ggraph`
- Saves and returns the following output:
  - Meta-data for the nodes in the network, including network properties
  - Meta-data for the clusters in the network
  - Network graph plot (returned as `ggraph` object and saved as pdf)
  - Network adjacency matrix in sparse matrix format
  - `igraph` object containing the list of edges in the network graph

See `vignette("buildRepSeqNetwork")` for a tutorial on the function’s
usage and output.

# Visualization

The `buildRepSeqNetwork` function includes various arguments that
facilitate customization of the network visualization. See
`vignette("network_visualization")` for an overview.

# Searching for Associated Clusters

Given multiple samples of AIRR-Seq data, the `NAIR` package can be used
to search for TCR/BCR clusters associated with a binary variable of
interest, such as a disease condition, treatment or clinical outcome.

See `vignette("associated_clusters")` for a detailed tutorial.

# Searching for Public Clusters

The `NAIR` package includes a set of functions that facilitate searching
for public TCR/BCR clusters across multiple samples of AIRR-seq data.

In this context, a public cluster consists of similar TCR/BCR clones
(e.g., those whose CDR3 amino acid sequences differ by at most one amino
acid) that are shared across samples (e.g., across individuals or across
time points for a single individual).

See `vignette("public_clusters")` for a detailed tutorial.
