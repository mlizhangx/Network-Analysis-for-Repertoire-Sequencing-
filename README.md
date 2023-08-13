
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-/actions/workflows/check-standard.yaml)
<!-- badges: end -->

# NAIR: Network Analysis of Immune Repertoire

We introduce the `NAIR` package, a powerful tool for analyzing the
adaptive immune repertoire using network analysis based on similarities
among receptor sequences, which implements methods from the following
paper:

[Hai Yang, Jason Cham, Brian Neal, Zenghua Fan, Tao He and Li Zhang.
(2023). NAIR: Network Analysis of Immune Repertoire. *Frontiers in
Immunology*, vol. 14.
https://doi.org/10.3389/fimmu.2023.1181825](https://www.frontiersin.org/articles/10.3389/fimmu.2023.1181825/full)

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

The development version of `NAIR` is [hosted on
Github](https://github.com/mlizhangx/Network-Analysis-for-Repertoire-Sequencing-)
and can be installed from source using the following command:

``` r
devtools::install_github(
  "mlizhangx/Network-Analysis-for-Repertoire-Sequencing-",
  dependencies = TRUE, 
  build_vignettes = TRUE
)
```

<!-- Installing the development version requires a toolchain compiler. On Windows, this means downloading and installing Rtools. On MacOS, this entails installing XCode Command Line Tools ("XCode CLI") and the correct version of gfortran for your macOS version (instructions [here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/)). -->

# Getting Started

[The website for the `NAIR`
package](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/)
contains articles to help users get started. All package documentation
and vignettes can be found there.

Once the package is installed, package vignettes and documentation can
be accessed offline from within R using the following commands:

``` r
# Package vignette
vignette("NAIR", package = "NAIR")

# Display vignettes in html browser
browseVignettes("NAIR")
```

## The `buildRepSeqNetwork` function

General network analysis on AIRR-Seq data is performed using
`buildRepSeqNetwork()`. This function does the following:

- Filters the input AIRR-Seq data according to user specifications
- Builds the network graph for the immune repertoire
- Performs additional network analysis, which can include:
  - Cluster analysis
  - Computation of network properties
- Generates customizable visual plots of the network graph
- Saves and returns the following output:
  - Metadata for the nodes in the network
  - Metadata for the clusters in the network
  - The plots of the network graph
  - The network graph itself, both as an adjacency matrix and as an
    `igraph` object

See [this
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/buildRepSeqNetwork.html)
for a tutorial on the function’s usage and output.

## Visualization

`buildRepSeqNetwork()` accepts various arguments for customizing the
network visualization. See [this
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/network_visualization.html)
for an overview.

## Searching for Associated Clusters

Given multiple samples of AIRR-Seq data, the `NAIR` package can be used
to search for TCR/BCR clusters associated with a binary variable of
interest, such as a disease condition, treatment or clinical outcome.

See [this
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/associated_clusters.html)
for a detailed tutorial.

## Searching for Public Clusters

The `NAIR` package includes a set of functions that facilitate searching
for public TCR/BCR clusters across multiple samples of AIRR-seq data.

In this context, a public cluster consists of similar TCR/BCR clones
(e.g., those whose CDR3 amino acid sequences differ by at most one amino
acid) that are shared across samples (e.g., across individuals or across
time points for a single individual).

See [this
vignette](https://mlizhangx.github.io/Network-Analysis-for-Repertoire-Sequencing-/articles/public_clusters.html)
for a detailed tutorial.
