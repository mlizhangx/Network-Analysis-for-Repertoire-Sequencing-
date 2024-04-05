
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) -->
<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
“[![metacran
downloads](https://cranlogs.r-pkg.org/badges/grand-total/NAIR)](https://cran.r-project.org/package=NAIR)”
[![R-CMD-check](https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-/actions/workflows/check-standard.yaml)
<!-- badges: end -->

# NAIR: Network Analysis of Immune Repertoire

`NAIR` is an R package for analyzing the adaptive immune repertoire
using network analysis based on similarities among receptor sequences.
It implements methods from the following paper:

[Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang.
(2023). NAIR: Network Analysis of Immune Repertoire. *Frontiers in
Immunology*, vol. 14.
https://doi.org/10.3389/fimmu.2023.1181825](https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full)

`NAIR` allows users to perform network analysis on Adaptive Immune
Receptor Repertoire Sequencing (AIRR-Seq) data, including computing
local and global network properties of nodes and clusters, which can
provide insights into the structural organization of the immune
repertoire network.

`NAIR` also enables users to search across multiple AIRR-Seq samples for
clones/clusters associated with subject characteristics, disease
conditions or clinical outcomes, as well as identify public
clones/clusters. This can help researchers identify potentially
important TCR/BCR clones.

To aid in interpretation of the immune repertoire network, `NAIR`
includes convenient functionality for generating customized network
visualizations.

#### What data does `NAIR` support?

`NAIR` supports bulk and single-cell immune repertoire sequence data for
T-cell or B-cell receptors (TCR or BCR).

- **Single-cell data:** Each row is a single cell
- **Bulk data:** Each row is a distinct TCR/BCR clone (unique
  combination of V-D-J genes and nucleotide sequence) and typically
  includes a corresponding measurement of clonal abundance (e.g., clone
  count and clone frequency/fraction)

#### How does `NAIR` model the immune repertoire as a network?

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
  colored according to desired metadata (e.g., disease status, sample,
  cluster, clonal abundance, etc.)

# Installation

To install the latest release version of `NAIR`, use the following
command:

``` r
install.packages("NAIR")
```

To install the latest development version of `NAIR` from source (which
requires compilation), use the following command:

``` r
devtools::install_github(
  "mlizhangx/Network-Analysis-for-Repertoire-Sequencing-",
  dependencies = TRUE, 
  build_vignettes = TRUE
)
```

<!-- Installing the development version requires a toolchain compiler. On Windows, this means downloading and installing Rtools. On MacOS, this entails installing XCode Command Line Tools ("XCode CLI") and the correct version of gfortran for your macOS version (instructions [here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/)). -->

# Getting Started

### Main Function

General network analysis on AIRR-Seq data is performed using
`buildRepSeqNetwork()` or its convenient alias `buildNet()`. This
function does the following:

- Filters the AIRR-Seq data according to user specifications
- Builds the network graph for the immune repertoire
- Performs additional network analysis, which can include:
  - Cluster analysis
  - Network properties
  - Customizable visual plots of the network graph
- Returns (and optionally saves) the following output:
  - The network graph (as `igraph` and adjacency matrix)
  - Metadata for the network
  - Metadata for the nodes in the network
  - Metadata for the clusters in the network
  - Plots of the network graph

See [this
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/buildRepSeqNetwork.html)
for a tutorial.

### Searching for Associated Clusters

Given multiple samples of bulk AIRR-Seq data, `NAIR` can be used to
search for TCR/BCR clusters associated with a binary variable of
interest, such as a disease condition, treatment or clinical outcome.
See [this
article](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html)
for a tutorial.

### Searching for Public Clusters

The `NAIR` package includes a set of functions that facilitate searching
for public TCR/BCR clusters across multiple samples of bulk AIRR-seq
data. In this context, a public cluster consists of similar TCR/BCR
clones (e.g., those whose CDR3 amino acid sequences differ by at most
one amino acid) that are shared across samples (e.g., across individuals
or across time points for a single individual). See [this
article](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html)
for a tutorial.

# Additional Resources

### Visualization

[This
article](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/network_visualization.html)
provides an introduction to the creation and customization of network
visualizations using `NAIR`.

### Network Properties and Cluster Analysis

[This
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/node_properties.html)
provides an introduction to computing node-level network properties with
`NAIR`.

[This
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/cluster_analysis.html)
explains how to perform cluster analysis with `NAIR`.

### Supplementary Functions

[This
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/supplementary.html)
provides an overview of `NAIR` utility functions that supplement the
main function `buildNet()`.
