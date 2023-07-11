
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
# vignettes
browseVignettes("NAIR")

# documentation files
help(package = "NAIR")
```

# The `buildRepSeqNetwork` function

General network analysis on AIRR-Seq data is performed using the
`buildRepSeqNetwork` function. This function does the following:

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

## Simulate Data for Demonstration

For demonstration purposes, we simulate some toy data using built-in
package functions.

We simulate a sample consisting of 200 observations (rows).

``` r
library(NAIR)
dir_out <- tempdir()
toy_data <- simulateToyData()
head(toy_data)
#>        CloneSeq CloneFrequency CloneCount SampleID
#> 1 TTGAGGAAATTCG    0.007873775       3095  Sample1
#> 2 GGAGATGAATCGG    0.007777102       3057  Sample1
#> 3 GTCGGGTAATTGG    0.009094910       3575  Sample1
#> 4 GCCGGGTAATTCG    0.010160859       3994  Sample1
#> 5 GAAAGAGAATTCG    0.009336593       3670  Sample1
#> 6 AGGTGGGAATTCG    0.010369470       4076  Sample1
```

``` r
nrow(toy_data)
#> [1] 200
```

## Primary Inputs

The `buildRepSeqNetwork` function has two required arguments, which
together specify the input data.

- The function’s first argument `data` accepts a data frame containing
  the AIRR-seq data, where each row corresponds to a single TCR/BCR
  clone (bulk data) or cell (single-cell data).
- The second argument `seq_col` specifies the data column containing the
  receptor sequences to be used as the basis of similarity between two
  cells or clones. This argument accepts either the column name (e.g.,
  `seq_col = "CloneSeq"`) or the column index (e.g., `seq_col = 1`,
  denoting the first column).

The `buildRepSeqNetwork` function can be executed using only these first
two arguments, as demonstrated below.

``` r
output <- buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", 
                             output_dir = dir_out)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Generating graph plot...
#>  Done.
#> Node-level meta-data saved to file:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/MyRepSeqNetwork_NodeMetadata.csv
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" style="display: block; margin: auto;" />

    #> Network graph plots saved to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/MyRepSeqNetwork.pdf
    #> Network igraph saved in edgelist format to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/MyRepSeqNetwork_EdgeList.txt
    #> Adjacency matrix saved to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/MyRepSeqNetwork_AdjacencyMatrix.mtx

A visual plot of the network graph is automatically printed to the
plotting window in R.

The function returns a list containing the following items:

``` r
names(output)
#> [1] "igraph"           "adjacency_matrix" "node_data"        "plots"
```

The item `node_data` is a data frame where each row corresponds to a
node in the network graph:

``` r
head(output$node_data)
#>         CloneSeq CloneFrequency CloneCount SampleID
#> 2  GGAGATGAATCGG    0.007777102       3057  Sample1
#> 5  GAAAGAGAATTCG    0.009336593       3670  Sample1
#> 8  GGGGAGAAATTGG    0.006220155       2445  Sample1
#> 11 GGGGGAGAATTGC    0.012969469       5098  Sample1
#> 12 GGGGGGGAATTGC    0.009079646       3569  Sample1
#> 13 AGGGGGAAATTGG    0.014941093       5873  Sample1
```

The output data frame includes all of the columns that were present in
the input data frame. In some cases, the user may not wish for all
columns to be included. In this case, the user can specify the columns
to be included using the `subset_cols` argument, which accepts a vector
containing either the column names or the column indices of the columns
to be kept. The column specified by the `seq_col` argument will
automatically be included regardless of the value of `subset_cols`.

## Input Filtering

The following options are useful for removing noise or irrelevant data
from the analysis, potentially improving the quality of the network
graph and downstream analysis results.

- `min_seq_length`: used to filter out TCR/BCR sequences by minimum
  length. This parameter sets a minimum sequence length, and any data
  rows with sequence lengths below this value will be dropped before
  constructing the network graph. By default, the minimum sequence
  length is set to 3.
- `drop_matches`: used to filter out TCR/BCR sequences by content. This
  parameter takes a character string or regular expression and checks
  each TCR/BCR sequence for a match. Data rows with matches are dropped
  before constructing the network graph.

## Settings When Building the Network

The settings used to construct the network can be customized using the
arguments below.

### Distance Function

The default method for measuring the similarity between TCR/BCR
sequences is the Hamming distance. It calculates the number of
differences between two sequences of the same length. If the sequences
have different lengths, the shorter sequence is extended by adding
non-matching characters to make it the same length as the longer
sequence.

The Levenshtein distance can be used as an alternative measurement to
determine the similarity between sequences. It calculates the minimum
number of single-character edits (insertions, deletions and
transformations) needed to transform one sequence into the other. This
method is particularly useful for comparing sequences of different
lengths and can account for insertions and deletions. When constructing
a network based on the similarity of CDR-3 nucleotide sequences, it is
preferable to use the Levenshtein distance instead of the default
Hamming distance by setting the argument `dist_type = "levenshtein"`.
However, the Levenshtein distance requires significantly more
computational time than the Hamming distance, which may be challenging
when working with large data sets having long TCR/BCR sequences.

### Distance Cutoff

The distance function specified in the `dist_type` argument (Hamming
distance by default) is used to quantify the similarity between TCR/BCR
sequences. The chosen distance measurement determines the distance
between two nodes in the network graph.

By default, two nodes in the graph are connected by an edge if their
distance is at most 1. However, if users want to adjust this cutoff, the
`dist_cutoff` argument can be set to a different value.

For example, if `dist_cutoff = 2`, then two nodes will be connected by
an edge if their distance is at most 2. The cutoff value controls the
stringency of the network construction and affects the number and
density of edges in the network.

### Isolated Nodes

By default, if a given node has no edges connecting it to any other
nodes in the network, it will be removed from the network graph and will
not be included in the output.

Notice that while our input data contained 200 rows, the `node_data`
data frame in the returned output contains fewer rows than this:

``` r
nrow(output$node_data)
#> [1] 122
```

These 122 rows correspond to the nodes that are joined by an edge to at
least one other node in the network.

If desired, all nodes can be kept in the network, including those that
do not have any edge connections to other nodes. This is accomplished by
setting the `drop_isolated_nodes` argument to `FALSE`.

## Network Properties and Cluster Analysis

The `buildRepSeqNetwork` function can perform additional network
analysis after constructing the network. This includes cluster analysis
and computation of network properties. Cluster analysis partitions the
network graph into densely-connected subgraphs, while network properties
describe the structural organization of the network.

### Node-Level Network Properties

Node-level network properties are properties that pertain to each
individual node in the network graph.

Some are local properties, meaning that their value for a given node
depends only on a subset of the nodes in the network. One example is the
network degree of a given node, which represents the number of other
nodes that are directly joined to the given node by an edge connection.

Other properties are global properties, meaning that their value for a
given node depends on all of the nodes in the network. An example is the
authority score of a node, which is computed using the entire graph
adjacency matrix (if we denote this matrix by $A$, then the principal
eigenvector of $A_T A$ represents the authority scores of the network
nodes).

The names of the node-level network properties that can be computed are
listed below:

- `degree`
- `cluster_id`
- `transitivity`
- `closeness`
- `centrality_by_closeness`
- `eigen_centrality`
- `centrality_by_eigen`
- `betweenness`
- `centrality_by_betweenness`
- `authority_score`
- `coreness`
- `page_rank`

The `cluster_id` property warrants further explanation and is discussed
in the [subsection on cluster analysis](#cluster-analysis).

#### Enabling Computation of Node-Level Properties

Use the setting `node_stats = TRUE` to enable computation of node-level
network properties.

The node data now contains node-level network properties in addition to
the original data columns:

``` r
names(output$node_data)
#>  [1] "CloneSeq"                  "CloneFrequency"           
#>  [3] "CloneCount"                "SampleID"                 
#>  [5] "degree"                    "transitivity"             
#>  [7] "eigen_centrality"          "centrality_by_eigen"      
#>  [9] "betweenness"               "centrality_by_betweenness"
#> [11] "authority_score"           "coreness"                 
#> [13] "page_rank"
```

``` r
head(output$node_data[ , c("CloneSeq", "degree", "authority_score")])
#>         CloneSeq degree authority_score
#> 2  GGAGATGAATCGG      1      0.00000000
#> 5  GAAAGAGAATTCG      3      0.00000000
#> 8  GGGGAGAAATTGG      2      0.04558649
#> 11 GGGGGAGAATTGC      4      0.15055366
#> 12 GGGGGGGAATTGC     10      0.52691798
#> 13 AGGGGGAAATTGG      5      0.14682343
```

By default, all of the available node-level properties are computed
except for `closeness`, `centrality_by_closeness` and `cluster_id`.

#### Specifying the Set of Node-Level Properties

The set of node-level properties that are computed can be specified
using the `stats_to_include` argument.

To compute all node-level properties, the user can simply specify
`stats_to_include = "all"`.

To specify a particular subset of the available node-level properties,
the user must provide a named logical vector following a particular
format. This vector can be created using the `chooseNodeStats` function.
Each node-level property name is an argument of `chooseNodeStats`, with
the argument accepting either `TRUE` or `FALSE` to specify whether the
property is computed. The default value for each argument agrees with
the default set of node properties seen in the previous example. In
other words, setting `stats_to_include = chooseNodeStats()` is the same
as leaving the `stats_to_include` argument unspecified.

Below is an example where the user wishes to compute the `closeness`
property in addition to the default properties, with the `page_rank`
property excluded:

``` r
# Modifying the default set of node-level properties
buildRepSeqNetwork(
  toy_data, "CloneSeq", node_stats = TRUE,
  stats_to_include = chooseNodeStats(closeness = TRUE, page_rank = FALSE))
```

If the user wishes to include only a small number of properties and
exclude the rest, this requires setting many argument values of
`chooseNodeStats` to `FALSE`, which can be inconvenient. In this case,
it may instead be easier to use the `exclusiveNodeStats` function, which
behaves in the same manner as `chooseNodeStats`, but all of its
arguments are set to `FALSE` by default.

``` r
# Include only the node-level properties specified below
buildRepSeqNetwork(
  toy_data, "CloneSeq", node_stats = TRUE, 
  stats_to_include = exclusiveNodeStats(degree = TRUE, transitivity = TRUE))
```

### Cluster Analysis

Cluster analysis involves using a community-finding algorithm to
partition the network graph into clusters (densely-connected subgraphs).

The cluster membership of each node is then recorded as the `cluster_id`
node-level network property, included as a variable in the node-level
metadata.

Cluster-level network properties of each cluster are computed. An
additional data frame containing these cluster-level network properties
will then be included in the output list returned by
`buildRepSeqNetwork`.

#### Enabling Cluster Analysis

To perform cluster analysis, use the setting `cluster_stats = TRUE`.
This will also compute the cluster-level network properties of each
cluster.

If the user has set `node_stats = TRUE` and has specified inclusion of
the `cluster_id` node-level property through the `stats_to_include`
argument, then cluster analysis will automatically be performed and the
`cluster_id` values will be added to the node-level data. However,
cluster-level network properties will not be computed unless the user
has set `cluster_stats = TRUE`.

#### Clustering Algorithm

By default, clustering is performed using the `cluster_fast_greedy`
algorithm from the `igraph` package. Other clustering algorithms from
the `igraph` package can be used instead of the default algorithm. The
algorithm is specified using the `cluster_fun` argument, which accepts
one of the following functions:

- `cluster_fast_greedy`
- `cluster_edge_betweenness`
- `cluster_fluid_communities`
- `cluster_infomap`
- `cluster_label_prop`
- `cluster_leading_eigen`
- `cluster_leiden`
- `cluster_louvain`
- `cluster_optimal`
- `cluster_spinglass`
- `cluster_walktrap`

For example, setting `cluster_fun = cluster_leiden` performs clustering
using the `cluster_leiden` algorithm.

For more information about a particular algorithm, users can refer to
its help documentation file. For example, the command
`?cluster_fast_greedy` loads the documentation file for the
`cluster_fast_greedy` (the default) algorithm, assuming the `NAIR`
package has been loaded (e.g., using `library(NAIR)`).

#### Cluster-Level Network Properties

If the user has set `cluster_stats = TRUE`, then after partitioning the
network into clusters, various network properties are computed for each
cluster.

These cluster-level network properties are recorded in their own data
frame, containing one row per cluster.

Below, we re-generate the output of `buildRepSeqNetwork`, this time with
cluster analysis included.

``` r
output <- buildRepSeqNetwork(toy_data, "CloneSeq", node_stats = TRUE, 
                             cluster_stats = TRUE, print_plots = FALSE,
                             output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Computing cluster membership within the network... Done.
#> Computing statistics for the 20 clusters in the network... Done.
#> Generating graph plot with nodes colored by cluster_id... Done.
```

The output list now includes an additional data frame containing the
cluster-level meta data:

``` r
names(output)
#> [1] "igraph"           "adjacency_matrix" "node_data"        "cluster_data"    
#> [5] "plots"
```

``` r
names(output$cluster_data)
#>  [1] "cluster_id"                  "node_count"                 
#>  [3] "mean_seq_length"             "mean_degree"                
#>  [5] "max_degree"                  "seq_w_max_degree"           
#>  [7] "agg_count"                   "max_count"                  
#>  [9] "seq_w_max_count"             "diameter_length"            
#> [11] "global_transitivity"         "assortativity"              
#> [13] "edge_density"                "degree_centrality_index"    
#> [15] "closeness_centrality_index"  "eigen_centrality_index"     
#> [17] "eigen_centrality_eigenvalue"
```

``` r
head(output$cluster_data[ , 1:6])
#>   cluster_id node_count mean_seq_length mean_degree max_degree seq_w_max_degree
#> 1          1         14           13.00        3.36          9    AAAAAAAAATTGC
#> 2          2         28           12.96        8.43         18    GGGGGGGAATTGG
#> 3          3          9           12.67        2.22          4     AGAAGAAAATTC
#> 4          4          6           13.00        3.33          9    GGGGGGAAATTGG
#> 5          5          6           12.00        2.17          3     AGGGAGGAATTC
#> 6          6         25           12.00        4.60         10     AAAAAAAAATTG
```

A brief description of each cluster-level property is given below:

- `node_count`: The number of nodes in the cluster.
- `mean_seq_length`: The mean sequence length in the cluster.
- `mean_degree`: The mean network degree in the cluster.
- `max_degree`: The maximum network degree in the cluster.
- `seq_w_max_degree`: The receptor sequence possessing the maximum
  degree within the cluster.
- `agg_count`: The aggregate count among all nodes in the cluster (based
  on the counts in `count_col`, if provided).
- `max_count`: The maximum count among all nodes in the cluster (based
  on the counts in `count_col`, if provided).
- `seq_w_max_count`: The receptor sequence possessing the maximum count
  within the cluster.
- `diameter_length`: The longest geodesic distance in the cluster.
- `assortativity`: The assortativity coefficient of the cluster’s graph,
  based on the degree (minus one) of each node in the cluster (with the
  degree computed based only upon the nodes within the cluster).
- `global_transitivity`: The transitivity (i.e., clustering coefficient)
  for the cluster’s graph, which estimates the probability that adjacent
  vertices are connected.
- `edge_density`: The number of edges in the cluster as a fraction of
  the maximum possible number of edges.
- `degree_centrality_index`: The cluster-level centrality index based on
  degree within the cluster graph.
- `closeness_centrality_index`: The cluster-level centrality index based
  on closeness, i.e., distance to other nodes in the cluster.
- `eigen_centrality_index`: The cluster-level centrality index based on
  the eigenvector centrality scores, i.e., values of the principal
  eigenvector of the adjacency matrix for the cluster.
- `eigen_centrality_eigenvalue`: The eigenvalue corresponding to the
  principal eigenvector of the adjacency matrix for the cluster.

#### Specifying the Count Column

Some cluster-level network properties, such as `agg_count` and
`max_count`, are only computed if the user specifies a column of the
input data containing counts for each row (i.e., clone count for bulk
data or Unique Molecular Identifier count for single-cell data). This
column is specified using the `count_col` function, which accepts a
column name or column index.

## Visualization

The `buildRepSeqNetwork` function includes various arguments that
facilitate customization of the network visualization.

### Node Colors

#### Color nodes according to meta data

The nodes in the graph can be colored according to a variable in the
node-level metadata. This is accomplished using the `color_nodes_by`
argument, which accepts a character string specifying the name of the
column containing the variable. The user can specify any variable that
will be present in the node-level meta data returned by
`buildRepSeqNetwork`, including any node-level network properties to be
computed.

For example, when calling `buildRepSeqNetwork`, we can color the nodes
based on the `CloneCount` column of the original input data by setting
`color_nodes_by = "CloneCount"` as seen below.

``` r
buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", node_stats = TRUE,
                   color_nodes_by = "CloneCount", output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Generating graph plot with nodes colored by CloneCount...
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

If the `color_nodes_by` argument is left unspecified,
`buildRepSeqNetwork` will attempt to use an available variable to color
the nodes. If a variable for clone count is provided using the
`count_col` argument, this will be used. Otherwise the `cluster_id`
variable will be used if it exists, followed by `degree`, then by other
network properties. If no suitable options are available, the nodes will
be left uncolored.

If the user does not wish the nodes to be colored dynamically, this can
be specified by setting `color_nodes_by = NULL`.

#### Adjust node color palette

A preset color palette can be used to color the nodes by providing the
appropriate character string to the `color_scheme` argument. The
argument accepts the following values:

- `"default"` for default `ggplot2` colors
- One of the following color scales from the `viridisLite` package.
  These scales are designed to maintain their perceptual uniformity
  (consistently across the scale, values close together appear similar,
  while values far apart appear distinct) when printed in grey scale or
  when viewed by individuals with common forms of color-blindness or
  color vision deficiency. More information, including images of the
  color scales, can be referenced in the [Introduction to the viridis
  color
  maps](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
  vignette hosted on CRAN’s website.
  - `"magma"` (or `"A"`)
  - `"inferno"` (or `"B"`)
  - `"plasma"` (or `"C"`)
  - `"viridis"` (or `"D"`)
  - `"cividis"` (or `"E"`)
  - `"rocket"` (or `"F"`)
  - `"mako"` (or `"G"`)
  - `"turbo"` (or `"H"`)
  - Any of the above color scales with `"-1"` appended, which reverses
    the direction of the color scale (e.g.,
    `color_scheme = "viridis-1"`)
- A discrete color palette from `grDevices::hcl.pals()` (these can only
  be used if the variable used to color the nodes is discrete)

Below we show an example of using the `plasma` color scheme from the
`viridisLite` package, with the direction of the color scale reversed:

``` r
# Using the "plasma" color scheme with reversed color gradient
buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", node_stats = TRUE,
                   color_nodes_by = "transitivity", 
                   color_scheme = "plasma-1", output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Generating graph plot with nodes colored by transitivity...
```

<img src="man/figures/README-unnamed-chunk-20-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

### Node Size

#### Adjust uniform node size

By default, all nodes are drawn at a uniform size value of 0.5. This
results in the nodes appearing very small when the plot is saved using
the default pdf dimensions of 12 inches wide by 10 inches tall. This
default behavior is intended to prevent nodes from overlapping and
obscuring edges in larger networks.

A different uniform node size can be specified by providing a numeric
value to the `size_nodes_by` argument. Below, we set the node size to
1.5, which is three times as large as the default node size:

``` r
buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", node_stats = TRUE,
                   color_nodes_by = "transitivity", color_scheme = "plasma-1",
                   size_nodes_by = 1.5, output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Generating graph plot with nodes colored by transitivity...
```

<img src="man/figures/README-unnamed-chunk-21-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

#### Size nodes according to meta data

Rather than draw all nodes using a uniform node size, it is possible to
for the nodes to be sized according to a variable within the node-level
metadata. This is achieved by providing the column name of the variable
to the `color_nodes_by` argument.

The minimum and maximum node sizes can be specified using the
`node_size_limits` argument, which accepts a numeric vector of length
two, where the first entry is the minimum node size and the second entry
is the maximum node size.

``` r
# Size nodes dynamically by network degree; use custom size range
buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", node_stats = TRUE,
                   color_nodes_by = "transitivity", color_scheme = "plasma-1",
                   size_nodes_by = "degree", 
                   node_size_limits = c(0.5, 1.5), output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Generating graph plot with nodes colored by transitivity...
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

### Labeling Clusters

If [cluster analysis is
performed](buildRepSeqNetwork.html#cluster-analysis) when calling
`buildRepSeqNetwork`, the output will include data on the resulting
clusters, with each cluster identified by a cluster ID.

In order to more easily reference these clusters within the visual plot
of the network graph, it is possible to label the clusters in the plot
with their cluster IDs. This must be done after calling
`buildRepSeqNetwork`, and is accomplished by using the the
`addClusterLabels` function to modify the plot contained in the output
of `buildRepSeqNetwork`. Note that `buildRepSeqNetwork` returns a list,
and one of the elements of this list is another list named `plots`,
which contains all plots generated by `buildRepSeqNetwork`.

The `addClusterLabels` function has two primary arguments. The plot to
be modified is provided to the `plot` argument. The entire output list
returned by `buildRepSeqNetwork` is provided to the `net` argument.

By default, only the 20 largest clusters by node count are labeled in
order to preserve legibility. This number can be changed by providing a
different value to the `top_n_clusters` argument. Instead of
prioritizing the clusters to label based on their node count, a
different variable within the cluster-level metadata can be used, so
long as the variable is numeric. This is achieved by providing the
column name of the variable to the `criterion` argument. The variable
should be present in the `cluster_data` data frame contained in the
output list of `buildRepSeqNetwork`. Rather than prioritizing the
clusters with the greatest values of this variable, those with the
lowest values can instead be prioritized for labeling by setting
`greatest_values = FALSE`.

The size of the cluster ID labels can be adjusted by providing a numeric
value to the `size` argument (the default is 5), and their color can be
adjusted by providing a valid character string to the `color` argument
(the default is `"black"`).

The `addClusterLabels` function assumes that the node-level meta data
includes a variable named `cluster_id` that records the cluster
membership of each node. This is the case when `buildRepSeqNetwork` is
called with cluster analysis enabled. If this variable has a name
different from `cluster_id` (e.g., if columns are manually renamed), the
correct variable name must be provided to the `cluster_id_col` argument.

``` r
# Generate network, but don't print the graph plot yet
network <- buildRepSeqNetwork(
  toy_data, seq_col = "CloneSeq", 
  node_stats = TRUE, cluster_stats = TRUE,
  color_nodes_by = "transitivity", color_scheme = "plasma-1",
  size_nodes_by = "degree", node_size_limits = c(0.5, 1.5), 
  print_plots = FALSE, output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Computing cluster membership within the network... Done.
#> Computing statistics for the 20 clusters in the network... Done.
#> Generating graph plot with nodes colored by transitivity... Done.

# Add labels to the two largest clusters and print the plot
addClusterLabels(plot = network$plots$transitivity,
                 net = network,
                 top_n_clusters = 2,
                 criterion = "node_count" # (the default)
                 )
#> Warning: Removed 120 rows containing missing values (`geom_text()`).
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" style="display: block; margin: auto;" />

### Generate Multiple Plots of a Network

The `buildRepSeqNetwork` can generate multiple plots of the same network
graph, with each plot coloring the nodes according to a different
variable. This is accomplished by providing a vector of column names to
the `color_nodes_by` argument instead of a single column name.

A different [color palette](#adjust-node-color-palette) can be used for
each plot. This is achieved by providing a character vector to the
`color_scheme` argument. This vector must have the same length as the
vector provided to the `color_nodes_by` argument. The color palette
specified by each entry will be applied to the corresponding plot. If,
instead, a single value is provided to the `color_scheme` argument, the
specified color palette will be applied to all of the plots.

``` r
buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", node_stats = TRUE, 
                   color_nodes_by = c("transitivity", "CloneCount"),
                   color_scheme = c("plasma-1", "default"),
                   size_nodes_by = "CloneCount", node_size_limits = c(0.1, 2.5),
                   output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Generating graph plot with nodes colored by transitivity...
```

<img src="man/figures/README-unnamed-chunk-24-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.
    #> Generating graph plot with nodes colored by CloneCount...

<img src="man/figures/README-unnamed-chunk-24-2.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

### Title and Subtitle

The `output_name` argument is used in the names of files saved by
`buildRepSeqNetwork`. By default, it is also used as the title for each
plot generated by `buildRepSeqNetwork`. The default value of
`output_name` is `"MyRepSeqNetwork"`.

By default, the subtitle of each plot contains information about the
settings used to construct the network, including the values of the
`dist_type` and `dist_cutoff` arguments.

A custom plot title and subtitle can be specified using the `plot_title`
and `plot_subtitle` arguments, respectively. Either element can be
omitted from the plot by supplying a `NULL` value to the corresponding
argument, as shown below.

``` r
buildRepSeqNetwork(toy_data, seq_col = "CloneSeq", node_stats = TRUE,
                   color_nodes_by = "transitivity", color_scheme = "plasma-1",
                   size_nodes_by = "degree", node_size_limits = c(0.5, 2.5),
                   plot_title = NULL, plot_subtitle = NULL, output_dir = NULL)
#> Input data contains 200 rows.
#> Removing sequences with length fewer than 3 characters... Done. 200 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 122 nodes (after removing isolated nodes).
#> Computing node-level network statistics... Done.
#> Generating graph plot with nodes colored by transitivity...
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

### Legend

#### Excluding Part or All of the Legend

By default, if the nodes are colored or sized dynamically according to a
variable, a legend will be included in the plot showing the color scale
and/or size scale.

The color scale can be manually excluded from the legend by setting
`color_legend = FALSE`. Similarly, setting `size_legend = FALSE` will
exclude the size scale from the legend.

An exception to the default behavior occurs when the variable used to
color the nodes is discrete with more than 20 distinct values. In this
case the color scale is automatically excluded from the legend in order
to prevent it from taking up excessive space in the plot. The name of
the variable used to color the nodes is then reported in the plot’s
subtitle for reference.

#### Legend Titles

By default, the title shown above the color scale in the legend is the
name of the variable used to color the nodes. A custom title can be
specified by providing a character string to the `color_title` argument.
Similarly, the default title for the size scale is the name of variable
used to size the nodes, and a custom title can be specified using the
`size_title` argument. Each of the `color_title` and `size_title`
arguments can also be set to `NULL` in order to omit the respective
title.

If a vector is provided to the `color_nodes_by` argument in order to
[generate multiple plots](#generate-multiple-plots-of-a-network), then
the `color_title` will accept a character vector of matching length,
where each entry is the title for the color legend in the corresponding
plot.

### Edge width

By default, the edges in the plot are drawn with a width of 0.1.

A different width can be used for the edges by setting the `edge_width`
argument to a different value.

### Excluding Plots

The following arguments can be used to control whether plots are
generated and/or displayed in R:

- Use `print_plots = FALSE` to prevent plots from being printed to the R
  plotting window. Plots will still be generated, included in the output
  and (assuming a [valid output directory](#output-directory)) saved to
  a pdf as usual.
- Use `plots = FALSE` to prevent plots from being generated entirely.

## Output Settings

### Output Directory

- By default, `buildRepSeqNetwork` saves its output to the current
  working directory. However, users can specify a directory using the
  `output_dir` argument.
- The specified output directory will be created if it does not already
  exist.
- Setting `output_dir` to `NULL` will prevent any output from being
  written to file.

### Output Type

- By default, each component of the list returned by
  `buildRepSeqNetwork` is saved as its own file, with the node-level and
  cluster-level meta data saved as csv files.
- For better compression and fewer files, users can specify
  `output_type = "rds"` or `output_type = "rda"`, which will save the
  entire output list to a single rds file or a single rda file,
  respectively.

### Plots

- Regardless of the value of `output_type`, a separate pdf file
  containing the graph plot(s) is created in `output_dir`. The
  dimensions (in inches) for the pdf can be adjusted using the
  `plot_width` and `plot_height` arguments, with the defaults being `12`
  and `10`, respectively. The pdf file is not created if no plots are
  generated (i.e., if the argument `plots` is set to `FALSE`).
- **Note:** the `ggraph` object for each plot is only saved if the user
  sets `output_type = "rds"` or `output_type = "rda"`. Using one of
  these settings is recommended used if the user wishes to modify any
  plots in the future. Note, however, that plots can always be
  re-generated from the node-level meta data using the
  `generateGraphPlot` function.

### Output File Name(s)

- By default, the name of each saved file starts with `MyRepSeqNetwork`.
  This prefix can be changed by providing a character string to the
  `output_name` argument.
- If `output_type` is set to `"rds"` or `"rda"`, then the name of the
  file will be the value of the `output_name` argument followed by the
  appropriate file extension (either `.rds` or `.rda`).
- The name of the pdf file for the graph plot will be the value of the
  `output_name` argument followed by the `.pdf` file extension.

# Utility Functions

The `NAIR` package contains various utility functions designed to
complement the main `buildRepSeqNetwork` function. These utility
functions streamline various aspects of working with the output of
`buildRepSeqNetwork`, providing greater flexibility in the
implementation of customized and/or downstream analyses. Some of the
roles played by these utility functions include the following:

- Modifying or augmenting plots produced by the `buildRepSeqNetwork`
  function
- Using the output of the `buildRepSeqNetwork` function to generate new
  plots
- Performing cluster analysis or computing network properties for a
  network constructed using the `buildRepSeqNetwork` function.

For more details:

``` r
vignette(topic = "Utility Functions", package = "NAIR")
```

# Finding Associated Clones

Given multiple samples of AIRR-Seq data, the `NAIR` package can be used
to search for TCR/BCR clusters associated with a binary variable of
interest, such as a disease condition, treatment or clinical outcome.

We first provide a brief conceptual overview, followed by a
demonstration in which we explain the process in greater detail.

## Overview of Process

1.  **Identify associated sequences.** Divide the subjects into two
    groups based on the two levels of the binary variable. Identify
    TCR/BCR sequences that exhibit a statistically significant
    difference in frequency between the two groups using Fisher’s exact
    test.
2.  **Identify clones with sequences similar to the associated
    sequences.** For each associated sequence, all sequences that fall
    within a certain distance (e.g., those that differ by at most one
    amino acid) comprise its neighborhood. From all samples, identify
    all clones whose sequences belong to this neighborhood.
3.  **Construct global network using identified clones and perform
    clustering.** Combine the clones from all neighborhoods into a
    single global network. Perform cluster analysis and assign
    membership to the global clusters. These clusters are considered as
    the associated clusters.

## Simulate Data for Demonstration

We simulate some toy data for demonstration.

As our binary variable of interest, we consider a single treatment
factor with two levels, labeled treatment and control.

In each of the two groups, we simulate 15 samples, each containing 30
observations. The generation probabilities of the possible sequences are
fixed within each group. In order to simulate the treatment effect, the
generation probabilities of certain sequences differ substantially
between the two groups.

Each sample is saved in a separate file using the .rds file format. The
files are named “`Sample1.rds`”, “`Sample2.rds`”, etc. The file path of
their directory is saved to the R environment variable
`dir_input_samples` for later reference.

``` r
# Use temp dir
data_dir <- tempdir()

# Directory to store input files
dir_input_samples <- file.path(data_dir, "input_samples")
dir.create(dir_input_samples, showWarnings = FALSE)

# Number of samples by control/treatment group
samples_c <- samples_t <- 15 
samples <- samples_c + samples_t
sample_size <- 30 # (seqs per sample)          

# sequences (first five are chosen to be associated with treatment)
base_seqs <- c("CASSGAYEQYF", "CSVDLGKGNNEQFF", 
               "CASSIEGQLSTDTQYF", 
               "CASSEEGQLSTDTQYF",
               "CASSPEGQLSTDTQYF",
               "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF") 

# relative generation probabilities by control/treatment group
pgen_c <- matrix(rep(c(rep(1, 5), rep(30, 3)), times = samples_c), 
                 nrow = samples_c, byrow = TRUE)
pgen_t <- matrix(rep(c(1, 1, rep(1/3, 3), rep(2, 3)), times = samples_t), 
                 nrow = samples_t, byrow = TRUE)
pgen <- rbind(pgen_c, pgen_t)

# Simulate the data
library(NAIR)
simulateToyData(    
  samples = samples, sample_size = sample_size,
  prefix_length = 1, prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = base_seqs, affix_probs = pgen, num_edits = 0,
  output_dir = dir_input_samples, no_return = TRUE)
#> [1] TRUE
```

The first few rows of the data for the first sample appear as follows:

``` r
# View first few rows of data for sample 1
head(readRDS(file.path(dir_input_samples, "Sample1.rds")))
#>        CloneSeq CloneFrequency CloneCount SampleID
#> 1 CASDAGVFQPQHF     0.02606559       2832  Sample1
#> 2 CASDAGVFQPQHF     0.03718396       4040  Sample1
#> 3 CASSHRGTDTQYF     0.03182726       3458  Sample1
#> 4 CASDAGVFQPQHF     0.04615781       5015  Sample1
#> 5 RASSLAGNTEAFF     0.06006498       6526  Sample1
#> 6 CASDAGVFQPQHF     0.03363123       3654  Sample1
```

## 1. Find Associated Sequences

The first step is to conduct a systematic search for associated
sequences within the provided samples using the `findAssociatedSeqs`
function. This search is a two-stage procedure. The unique receptor
sequences are first filtered according to basic criteria in order to
narrow the list of candidates. Then for each candidate sequence, we
compute the P-value for Fisher’s exact test of independence between the
binary variable of interest and the observed presence of the sequence
within a sample/subject. The user specifies the P-value cutoff below
which an association is detected.

Below, we explain the usage and behavior of the `findAssociatedSeqs`
function.

### 1.1 Filter and Cutoff Settings

The `findAssociatedSeqs` function has several parameters that control
the filter criteria used to determine which sequences are considered for
testing, as well as an argument to specify the P-value cutoff below
which an association is detected. These arguments are presented below.

#### 1.1.1 Sample Membership

By default, only sequences that appear in at least 5 samples will be
considered. This can be changed by setting the `min_sample_membership`
argument to a different value. Setting the value to `NULL` bypasses this
check.

#### 1.1.2 Sequence Length

By default, only sequences that contain at least 7 characters will be
considered. This can be changed by setting the `min_seq_length` argument
to a different value. Setting the value to `NULL` bypasses this check.

#### 1.1.3 Sequence Content

Sequences containing characters `*` or `_` will be excluded from
consideration. This can be changed using the `drop_matches` argument,
which accommodates a character string or regular expression specifying
the pattern of content to search for. The content of each sequence is
checked for a match to this pattern using the `grep` function from base
R. If a match is found, the sequence is excluded from consideration.
Setting the value to `NULL` bypasses this check.

For details on how the pattern matching is performed, please refer to
the base R documentation files for `regex` and `grep`.

#### 1.1.4 P-value Cutoff

By default, sequences with a P-value below 0.05 on Fisher’s exact test
are included in the output of `findAssociatedSeqs`. However, users may
wish to impose a stronger burden of evidence to control the false
discovery rate under multiple testing. The cutoff can be set to a
different value using the `pval_cutoff` argument. The lower the cutoff
value, the stronger the evidence of an association is required for a
sequence to be included in the output.

It should be noted, however, that the sequences returned by
`findAssociatedSeqs` are ordered by P-value, and any subset of them can
be used in the following step (2). Thus, imposing a stricter P-value
cutoff can also be done indirectly after calling `findAssociatedSeqs`,
by subsetting the results according to P-value.

### 1.2 Input File List

The main argument of the `findAssociatedSeqs` function is the
`file_list` argument, which accepts a vector containing file paths. Each
path corresponds to a distinct AIRR-Seq data file representing an
individual sample.

Below, we prepare the vector `input_files` to be provided to the
`file_list` argument of `findAssociatedSeqs`:

``` r
# input files for step 1 (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))
head(input_files)
#> [1] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample1.rds"
#> [2] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample2.rds"
#> [3] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample3.rds"
#> [4] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample4.rds"
#> [5] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample5.rds"
#> [6] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample6.rds"
```

### 1.3 Input Type

The file format of the input files for `findAssociatedSeqs` is specified
using the `input_type` parameter. The supported formats include `"rds"`,
`"rda"`, `"csv"`, as well as files that can be read using the
`read.table` function, such as `"tsv"` and `"txt"`.

For text formats such as `"csv"`, `"tsv"` and `"txt"`, users can specify
the separation option by utilizing the `sep` argument. The default
setting `sep = ""` accommodates all forms of white space, i.e., one or
more spaces, tabs, newlines or carriage returns. In addition, it is
important to note that the first line of the data is assumed to be the
header by default. To disable this behavior and treat the first line as
data, users must set the `header` parameter to `FALSE`.

Our samples are stored in .rds files, so we use `input_type = "rds"`.

### 1.4 Specifying the Sequence Column

The `seq_col` argument is used to specify the column containing the
clone sequences in the input data for each sample. The argument accepts
either the column name or column index.

In our simulated data, the column containing the clone sequences is
named `CloneSeq`.

### 1.5 Assigning Samples to Groups

The `group_ids` argument is used to assign each sample to one of the two
groups, representing the two levels of the binary variable of interest.
The argument accepts a vector of the same length as `file_list`. Each
entry of `group_ids` is assigned as a group label to the sample in the
corresponding entry of `file_list`. Any values may be used for the group
labels, but the vector must contain exactly two unique values.

For instance, in our simulated data, the first half of the samples
belong to the control group, while the second half belong to the
treatment group. Thus, we should assign one group label to the first 15
samples and a different group label to the last 15 samples. Here we
choose to label the first 15 samples as `"reference"` and the last 15
samples as `"comparison"`. However, it should be noted that the results
will be unchanged if the labels are reversed, or if a different pair of
labels is used, as long as the first 15 samples are assigned to one
group and the last 15 samples to the other.

The vector we will provide to the `group_ids` argument of
`findAssociatedSeqs` is created below:

``` r
# group label assignments for the samples
group_labels <- c(rep("reference", samples_c), rep("comparison", samples_t))
group_labels
#>  [1] "reference"  "reference"  "reference"  "reference"  "reference" 
#>  [6] "reference"  "reference"  "reference"  "reference"  "reference" 
#> [11] "reference"  "reference"  "reference"  "reference"  "reference" 
#> [16] "comparison" "comparison" "comparison" "comparison" "comparison"
#> [21] "comparison" "comparison" "comparison" "comparison" "comparison"
#> [26] "comparison" "comparison" "comparison" "comparison" "comparison"
```

### 1.6 Assigning Samples to Subjects (If Applicable)

The `subject_ids` argument can be used to assign each sample to a
particular subject. The argument accepts a vector of the same length as
`file_list`. Each entry of `subject_ids` is assigned as a subject ID to
the sample in the corresponding entry of `file_list`.

If the `subject_ids` argument is omitted, Fisher’s exact test treats
each sample as an independent observational unit. In this case, the
relevant contingency table involves counts of **samples** possessing a
given sequence.

If subject IDs are provided, each subject’s collection of samples is
treated as a single observational unit. The relevant contingency table
then involves counts of **subjects** possessing a given sequence. This
allows a sequence to be counted at most once per subject, and results in
each subject being counted exactly once in each margin.

Subject IDs should be provided when the binary variable of interest is
subject-specific and the data contains multiple samples from a single
subject. Subject IDs should be omitted when the binary variable of
interest is sample-specific or each sample comes from a different
subject.

### 1.7 Specifying the Clone Frequency Column (Optional)

The `freq_col` argument can be used to specify a column containing clone
frequencies in the input data for each sample. The argument accepts
either the column name or column index.

If clone frequencies are provided, the maximum clone frequency (across
all samples) for each associated sequence will be included in the
content of the data frame returned by `findAssociatedSeqs`.

### 1.8 Output Settings

The `findAssociatedSeqs` function returns a data frame containing the
associated sequences along with some additional information. The format
and contents of this data frame will be explained below after executing
the function.

By default, the data frame returned by `findAssociatedSeqs` is also
saved to the current working directory as a csv file named
`associated_seqs.csv`.

A different file name and/or directory can be specified by providing a
file path to the `outfile` argument. For example, setting
`outfile = "myfile.csv"` will save the file to the current working
directory as `myfile.csv`, while setting
`outfile = "~/myfolder/myfile.csv"` will save the file within the
subdirectory `myfolder` located within the current working directory.

The user can also specify `outfile = NULL` in order to prevent the
output from being saved.

### 1.9 Execution and Output

We execute the `findAssociatedSeqs` function using the inputs we
prepared earlier for the `file_list` and `group_ids` arguments:

``` r
# search across samples for associated sequences using Fisher's exact test
associated_seqs <- findAssociatedSeqs(
  file_list = input_files, input_type = "rds", 
  group_ids = group_labels, 
  seq_col = "CloneSeq", 
  min_seq_length = NULL, drop_matches = NULL, 
  min_sample_membership = NULL, 
  pval_cutoff = 0.1,
  outfile = NULL)
#> Data contains 30 samples, 15 of which belong to group reference and 15 of which belong to group comparison.
#> >>> Loading and compiling data from all samples:
#> Loading sample 1: Input data contains 30 rows.
#> Loading sample 2: Input data contains 30 rows.
#> Loading sample 3: Input data contains 30 rows.
#> Loading sample 4: Input data contains 30 rows.
#> Loading sample 5: Input data contains 30 rows.
#> Loading sample 6: Input data contains 30 rows.
#> Loading sample 7: Input data contains 30 rows.
#> Loading sample 8: Input data contains 30 rows.
#> Loading sample 9: Input data contains 30 rows.
#> Loading sample 10: Input data contains 30 rows.
#> Loading sample 11: Input data contains 30 rows.
#> Loading sample 12: Input data contains 30 rows.
#> Loading sample 13: Input data contains 30 rows.
#> Loading sample 14: Input data contains 30 rows.
#> Loading sample 15: Input data contains 30 rows.
#> Loading sample 16: Input data contains 30 rows.
#> Loading sample 17: Input data contains 30 rows.
#> Loading sample 18: Input data contains 30 rows.
#> Loading sample 19: Input data contains 30 rows.
#> Loading sample 20: Input data contains 30 rows.
#> Loading sample 21: Input data contains 30 rows.
#> Loading sample 22: Input data contains 30 rows.
#> Loading sample 23: Input data contains 30 rows.
#> Loading sample 24: Input data contains 30 rows.
#> Loading sample 25: Input data contains 30 rows.
#> Loading sample 26: Input data contains 30 rows.
#> Loading sample 27: Input data contains 30 rows.
#> Loading sample 28: Input data contains 30 rows.
#> Loading sample 29: Input data contains 30 rows.
#> Loading sample 30: Input data contains 30 rows.
#> All samples loaded.
#> Extracting list of unique sequences... Done.  8 unique sequences present.
#> Computing sample membership (this may take a while)... Done.
#> Filtering by Fisher's exact test P-value... Done. 4 sequences remain.
#> All done. Sorting results by Fisher's exact test P-value and returning.
```

`findAssociatedSeqs` returns a data frame containing the receptor
sequences found to be associated with the binary variable based on
Fisher’s exact test using the specified P-value cutoff. Each row
corresponds to a unique sequence and includes the following variables:

- `ReceptorSeq`: The unique receptor sequence
- `fisher_pvalue`: The P-value on Fisher’s exact test for independence
  between the receptor sequence and the binary variable of interest
- `shared_by_n_samples`: The number of samples in which the sequence was
  observed
- `samples_g0`: Of the samples in which the sequence was observed, the
  number of samples belonging to the first group (first unique value of
  `group_ids`)
- `samples_g1`: Of the samples in which the sequence was observed, the
  number of samples belonging to the second group (second unique value
  of `group_ids`)
- `shared_by_n_subjects`: The number of subjects in which the sequence
  was observed (only present if subject IDs are specified through
  `subject_ids`)
- `subjects_g0`: Of the subjects in which the sequence was observed, the
  number of subjects belonging to the first group (only present if
  subject IDs are specified through `subject_ids`)
- `subjects_g1`: Of the subjects in which the sequence was observed, the
  number of subjects belonging to the second group (only present if
  subject IDs are specified through `subject_ids`)
- `label`: A character string summarizing the above information. Also
  includes the maximum in-sample clone frequency across all samples, if
  available.

``` r
# view first few rows of output
head(associated_seqs)
#>        ReceptorSeq fisher_pvalue shared_by_n_samples samples_g0 samples_g1
#> 8   CSVDLGKGNNEQFF  1.052106e-05                  18          3         15
#> 7      CASSGAYEQYF  1.157316e-04                  17          3         14
#> 4 CASSEEGQLSTDTQYF  5.197401e-03                  10          1          9
#> 5 CASSIEGQLSTDTQYF  6.559548e-02                  16          5         11
#>                                                                                                                  label
#> 8 Sequence present in 18 samples (3 in group reference, 15 in group comparison)\nFisher's exact test P-value: 1.05e-05
#> 7 Sequence present in 17 samples (3 in group reference, 14 in group comparison)\nFisher's exact test P-value: 0.000116
#> 4    Sequence present in 10 samples (1 in group reference, 9 in group comparison)\nFisher's exact test P-value: 0.0052
#> 5   Sequence present in 16 samples (5 in group reference, 11 in group comparison)\nFisher's exact test P-value: 0.0656
```

The rows of the data frame are ordered by Fisher’s exact test $P$-value.

## 2. Find Associated Clones

In the [previous step (1)](#1.-find-associated-sequences), we used
`findAssociatedSeqs` to identify receptor sequences associated with the
binary variable of interest.

The next step is to use the `findAssociatedClones` function to search
across samples and identify all clones with sequences similar to the
associated sequences identified in step 1.

For each associated sequence, we define its neighborhood as the set of
all sequences that fall within a specified distance (e.g., a maximum
Hamming distance of 1). We then identify all clones (from all samples)
whose sequences belong to this neighborhood.

The data for each associated sequence’s neighborhood is then saved to a
separate file to be used as an input in [step
3](#3.-global-network-of-associated-clusters).

### 2.1 Specifying the Sample Data

In order for `findAssociatedClones` to conduct its search, we must
provide specifications for the sample data, just as we did when calling
`findAssociatedSeqs` in step 1. This is done using the arguments
`file_list`, `input_type`, `group_ids` and `seq_col`, which behave in
the same manner as seen earlier in the `findAssociatedSeqs` function.

### 2.2 Assigning Subject IDs and Sample IDs (Optional)

The `subject_ids` argument allows for subject IDs to be assigned to the
samples in the same manner as in the `findAssociatedSeqs` function. If
subject IDs are provided, each clone’s subject ID will be included in
the data for each associated sequence’s neighborhood.

The `sample_ids` argument allows for custom sample IDs to be assigned.
By default, the samples are labeled numerically according to the order
they appear in `file_list`. Each clone’s sample ID is included in the
data for each associated sequence’s neighborhood.

The `subject_ids` and `sample_ids` arguments both accept a vector of the
same length as `file_list`. Each entry of `subject_ids` (respectively,
`sample_ids`) is assigned as a subject ID (respectively, sample ID) to
the sample in the corresponding entry of `file_list`.

### 2.3 Specifying the Associated Sequences

The associated sequences are specified via the `assoc_seqs` argument,
which accepts a character vector.

Typically, the vector provided to `assoc_seqs` will be the `ReceptorSeq`
column of the data frame returned by `findAssociatedSeqs`. This
considers all of the associated sequences found in step 1.

However, it may be desirable to consider only a subset of the sequences
returned by `findAssociatedSeqs`. The sequences are ordered by Fisher’s
exact test P-value to facilitate reference. For example, if we had many
associated sequences, we could choose to consider only the 10 with the
lowest P-values by specifying
`assoc_seqs = associated_seqs$ReceptorSeq[1:10]`.

### 2.4 Neighborhood Distance Settings

By default, each associated sequence’s neighborhood includes all
sequences with a Hamming distance of at most 1 from the associated
sequence.

The type of distance metric and the distance threshold used to determine
the neighborhoods can be adjusted using the `dist_type` and `nbd_radius`
arguments.

For example, setting `dist_type = "lev"` and `nbd_radius = 2` results in
each neighborhood containing all sequences with a Levenshtein distance
of at most 2 from the associated sequence.

### 2.5 Output Settings

The `findAssociatedClones` function does not return any direct output.
Instead, it saves the network data for the associated sequence’s
neighborhoods to files that will be used as inputs in [step
3](#3.-global-network-of-associated-clusters).

The file path for the output directory is specified using the
`output_dir` argument. The output directory will be created if it does
not already exist.

One file is saved for each associated sequence. By default, each file is
saved as a csv file, but this can be changed using the `output_type`
argument. Other valid options include `"tsv"`, `"rds"` and `"rda"`.

``` r
# output directory for current step
dir_nbds <- file.path(data_dir, "assoc_seq_nbds")
```

### 2.6 Execution and Output

We execute the `findAssociatedClones` function using the inputs we
prepared earlier for the `file_list` and `group_ids` arguments:

``` r
# Identify clones in a neighborhood around each associated sequence
findAssociatedClones(
  file_list = input_files, input_type = "rds", 
  group_ids = group_labels, 
  seq_col = "CloneSeq", 
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL, drop_matches = NULL,
  output_dir = dir_nbds)
#> <<< Beginning search for associated clones >>>
#> Processing sample 1 of 30 (1):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 2 of 30 (2):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 3 of 30 (3):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 4 of 30 (4):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 5 of 30 (5):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 6 of 30 (6):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 7 of 30 (7):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 8 of 30 (8):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 9 of 30 (9):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 10 of 30 (10):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 11 of 30 (11):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 12 of 30 (12):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 13 of 30 (13):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 14 of 30 (14):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 15 of 30 (15):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 16 of 30 (16):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 17 of 30 (17):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 18 of 30 (18):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 19 of 30 (19):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 20 of 30 (20):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 21 of 30 (21):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 22 of 30 (22):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 23 of 30 (23):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 24 of 30 (24):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 25 of 30 (25):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 26 of 30 (26):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 27 of 30 (27):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 28 of 30 (28):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 29 of 30 (29):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> Processing sample 30 of 30 (30):
#> Input data contains 30 rows.
#> Finding clones in a neighborhood of each associated sequence... Done.
#> >>> Done processing samples. Compiling results:
#> Gathering data from all samples for sequence 1 (CSVDLGKGNNEQFF)... Done.
#> Gathering data from all samples for sequence 2 (CASSGAYEQYF)... Done.
#> Gathering data from all samples for sequence 3 (CASSEEGQLSTDTQYF)... Done.
#> Gathering data from all samples for sequence 4 (CASSIEGQLSTDTQYF)... Done.
#> >>> All tasks complete. Output is contained in the following directory:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/assoc_seq_nbds
```

The following files are created in the directory specified by the
`output_dir` argument:

``` r
# Files created by findAssociatedClones
list.files(dir_nbds)
#> [1] "CASSEEGQLSTDTQYF.csv" "CASSGAYEQYF.csv"      "CASSIEGQLSTDTQYF.csv"
#> [4] "CSVDLGKGNNEQFF.csv"
```

Each file contains the neighborhood data for a single associated
sequence.

## 3. Global Network of Associated Clusters

Now that we have identified the clones in each associated sequence’s
neighborhood, the final step is to use the
`buildAssociatedClusterNetwork` function to combine the clones from all
neighborhoods into a single global network. We then use clustering
analysis to partition the global network into clusters, which are
considered as the associated clusters.

### 3.1 Specifying the Neighborhood Data Files

The files created by `findAssociatedClones` in the [previous step
(2)](#2.-find-associated-clones) contain the data for each neighborhood.
These files are provided to `buildAssociatedClusterNetwork` by supplying
a character vector of file paths to the `file_list` argument. We create
this vector below.

``` r
# Files created by findAssociatedClones
nbd_files <- list.files(dir_nbds, full.names = TRUE)
```

### 3.2 Customization of Network Analysis

`buildAssociatedClusterNetwork` uses the same arguments as
[`buildRepSeqNetwork`](buildRepSeqNetwork.html) for customizing the
global network analysis.

#### 3.2.1 Distance Function and Distance Cutoff

By default, two nodes within the global network are joined by an edge if
the Hamming distance between their sequences is at most 1.

The Levenshtein distance can be used instead of the Hamming distance by
setting `dist_type = "lev"`.

The maximum distance for two nodes to be joined by an edge can be
changed using the `dist_cutoff` argument.

#### 3.2.2 Clustering Algorithm

After constructing the global network, `buildAssociatedClusterNetwork`
performs cluster analysis on the network nodes, partitioning the global
network graph into densely-connected subgraphs.

If the neighborhoods for two different associated sequences are very
similar, i.e., the network contains many edge connections joining the
nodes from one neighborhood to the nodes from the other neighborhood,
then the two neighborhoods will have a high chance of belonging to the
same cluster. Thus, cluster analysis assists in identifying distinct
groups of clones/sequences that are associated with the binary variable
of interest.

By default, clustering is performed using the `cluster_fast_greedy`
algorithm from the `igraph` package. Other clustering algorithms from
the `igraph` package can be used instead of the default algorithm. The
algorithm is specified using the `cluster_fun` argument, which accepts
one of the following functions:

- `cluster_fast_greedy`
- `cluster_edge_betweenness`
- `cluster_fluid_communities`
- `cluster_infomap`
- `cluster_label_prop`
- `cluster_leading_eigen`
- `cluster_leiden`
- `cluster_louvain`
- `cluster_optimal`
- `cluster_spinglass`
- `cluster_walktrap`

For example, setting `cluster_fun = cluster_leiden` performs clustering
using the `cluster_leiden` algorithm.

For more information about a particular algorithm, users can refer to
its help documentation file. For example, the command
`?cluster_fast_greedy` loads the documentation file for the
`cluster_fast_greedy` (the default) algorithm, assuming the `NAIR`
package has been loaded (e.g., using `library(NAIR)`).

### 3.3 Customization of Visual Plot

By default, the network graph plot produced by
`buildAssociatedClusterNetwork` colors the nodes according to the binary
variable of interest. This assists the user in validating each cluster’s
association to the binary variable. It also allows one to visually
distinguish clusters in which only a single level of the binary variable
is present (e.g., disease-only clusters).

If desired, a different variable can be used to color the nodes. This is
done by specifying the column for the variable to the `color_nodes_by`
argument, which accepts a column name or column index.

The `color_nodes_by` argument also accommodates a vector of column names
or a vector of column indices, in which case one plot will be created
for each column specified, with each plot coloring the nodes according
to the variable in its respective column. For example, setting
`color_nodes_by = c(3, 5)` will produce two plots per sample, where one
plot has the nodes colored according to the value in the third column of
the input data

Other arguments accepted by `buildRepSeqNetwork` to customize the
visualization can also be used when calling
`buildAssociatedClusterNetwork`.

### 3.4 Output Settings

The output returned by `buildAssociatedClusterNetwork` follows the same
format as the output of [`buildRepSeqNetwork`](buildRepSeqNetwork.html).
The function returns a list containing the node-level and cluster-level
meta data for the global network, as well as any plots generated, in
addition to the network adjacency matrix and the `igraph` network edge
list.

By default, the contents of the list returned by
`buildAssociatedClusterNetwork` are saved to the current working
directory. Each list element is saved as an individual file. The file
formats are the same default file formats used by
[`buildRepSeqNetwork`](buildRepSeqNetwork.html). In particular, the
node-level and cluster-level meta data are saved as csv files.

Alternatively, the user can save the entire output list to a single
compressed rds or rda file by setting `output_type = "rds"` or
`output_type = "rda"`, respectively.

By default, all files saved share the common file name prefix
`"AssociatedClusterNetwork"`. This common file name prefix can be set to
a different value by supplying a character string to the `output_name`
argument.

The output can be saved to a different directory by providing a file
path to the `output_dir` argument.

The user can also specify `output_dir = NULL` in order to prevent the
output from being saved.

### 3.5 Execution and Output

We execute the `buildAssociatedClusterNetwork` function using the input
we prepared earlier for the `file_list` argument:

``` r
# Combine neighborhoods and perform network analysis
all_clusters <- buildAssociatedClusterNetwork(
  file_list = nbd_files, 
  seq_col = "CloneSeq", 
  size_nodes_by = 1.5,
  output_dir = file.path(data_dir, "assoc_clusters"))
#> <<< Building network of associated clones >>>
#> Input data contains 202 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 202 nodes.
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 3 clusters in the network... Done.
#> Generating graph plot with nodes colored by GroupID...
#>  Done.
#> Node-level meta-data saved to file:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/assoc_clusters/AssociatedClusterNetwork_NodeMetadata.csv
#> Cluster-level meta-data saved to file:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/assoc_clusters/AssociatedClusterNetwork_ClusterMetadata.csv
```

<img src="man/figures/README-unnamed-chunk-37-1.png" width="100%" style="display: block; margin: auto;" />

    #> Network graph plots saved to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/assoc_clusters/AssociatedClusterNetwork.pdf
    #> Network igraph saved in edgelist format to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/assoc_clusters/AssociatedClusterNetwork_EdgeList.txt
    #> Adjacency matrix saved to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/assoc_clusters/AssociatedClusterNetwork_AdjacencyMatrix.mtx

In the visual plot of the network, we see that the treatment group shows
high representation in the clusters, affirming the association between
these sequences and the treatment variable.

The elements of the output list returned by
`buildAssociatedClusterNetwork` have the following names:

``` r
# output returned by buildAssociatedClusterNetwork
names(all_clusters)
#> [1] "igraph"           "adjacency_matrix" "node_data"        "cluster_data"    
#> [5] "plots"
```

The default files saved by `buildAssociatedClusterNetwork` have the
following names:

``` r
# Files saved by buildAssociatedClusterNetwork
list.files(file.path(data_dir, "assoc_clusters"))
#> [1] "AssociatedClusterNetwork.pdf"                
#> [2] "AssociatedClusterNetwork_AdjacencyMatrix.mtx"
#> [3] "AssociatedClusterNetwork_ClusterMetadata.csv"
#> [4] "AssociatedClusterNetwork_EdgeList.txt"       
#> [5] "AssociatedClusterNetwork_NodeMetadata.csv"
```

#### 3.5.1 Node-Level Meta Data

The `node_data` data frame contained in the output list contains the
following variables:

``` r
# variables in the node-level meta data
names(all_clusters$node_data)
#>  [1] "CloneSeq"                  "CloneFrequency"           
#>  [3] "CloneCount"                "SampleID"                 
#>  [5] "GroupID"                   "AssocSeq"                 
#>  [7] "degree"                    "cluster_id"               
#>  [9] "transitivity"              "eigen_centrality"         
#> [11] "centrality_by_eigen"       "betweenness"              
#> [13] "centrality_by_betweenness" "authority_score"          
#> [15] "coreness"                  "page_rank"
```

Notice that by default, all variables that were present in each sample’s
original data, such as `"CloneFrequency"` and `"CloneCount"`, are
automatically carried over into this data.

In addition, variables containing various node-level network properties
are also present.

#### 3.5.2 Cluster-Level Meta Data

The `cluster_data` data frame contained in the output list contains the
following variables:

``` r
# variables in the node-level meta data
names(all_clusters$cluster_data)
#>  [1] "cluster_id"                  "node_count"                 
#>  [3] "mean_seq_length"             "mean_degree"                
#>  [5] "max_degree"                  "seq_w_max_degree"           
#>  [7] "agg_count"                   "max_count"                  
#>  [9] "seq_w_max_count"             "diameter_length"            
#> [11] "global_transitivity"         "assortativity"              
#> [13] "edge_density"                "degree_centrality_index"    
#> [15] "closeness_centrality_index"  "eigen_centrality_index"     
#> [17] "eigen_centrality_eigenvalue"
```

Each row corresponds to a cluster in the global network, and each
variable corresponds to a cluster-level property.

### 3.6 Labeling the Global Clusters

In order to more easily cross-reference the clusters in the visual plot
with the clusters in the data, we can [label the clusters with their ID
numbers](network_visualization.html#labeling-clusters) as follows:

``` r
# Modify plot to add labels to the clusters
all_clusters$plots[[1]] <- 
  addClusterLabels(
    plot = all_clusters$plots[[1]],
    net = all_clusters,
    top_n_clusters = 3,
    criterion = "node_count", # (the default)
    size = 10
  )

# View modified plot
all_clusters$plots[[1]]
#> Warning: Removed 199 rows containing missing values (`geom_text()`).
```

<img src="man/figures/README-unnamed-chunk-42-1.png" width="100%" style="display: block; margin: auto;" />

### 3.7 Focusing on Individual Clusters of Interest

If we wish to focus on a particular cluster of interest within the
global network, we can build a network exclusively using the clones from
that cluster.

Below, we focus on the first cluster, which is also the largest cluster
by node count.

We generate two plots. In the first plot, we color each node according
to the receptor sequence of its corresponding clone. This provides a
more detailed account of the sequences that appear in the cluster and
their relative representation.

In the second plot, we color each node according to the sample in which
the corresponding clone originally appeared. This allows one to
distinguish whether the clones in the cluster come from many samples as
opposed to relatively few samples.

``` r
# focus on the first cluster
buildRepSeqNetwork(
  data = all_clusters$node_data[all_clusters$node_data$cluster_id == 1, ],
  seq_col = "CloneSeq", 
  color_nodes_by = c("CloneSeq", "SampleID"), 
  color_scheme = c("plasma", "turbo"),
  size_nodes_by = 3, 
  output_dir = NULL, output_name = "Cluster 1")
#> Input data contains 83 rows.
#> Removing sequences with length fewer than 3 characters... Done. 83 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 83 nodes (after removing isolated nodes).
#> Generating graph plot with nodes colored by CloneSeq...
```

<img src="man/figures/README-unnamed-chunk-43-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.
    #> Generating graph plot with nodes colored by SampleID...

<img src="man/figures/README-unnamed-chunk-43-2.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

# Finding Public Clones

The `NAIR` package includes a set of functions that facilitate searching
for public TCR/BCR clusters across multiple samples of AIRR-seq data.

In this context, a public cluster consists of similar TCR/BCR clones
(e.g., those whose CDR3 amino acid sequences differ by at most one amino
acid) that are shared across samples (e.g., across individuals or across
time points for a single individual).

We first provide a brief conceptual overview, followed by a
demonstration in which we explain the process in greater detail.

## Overview of Process

1.  **Filter the clusters within each sample’s network**. For each
    sample, construct the repertoire network and use cluster analysis to
    partition the nodes into clusters. Filter the data, keeping only
    those clusters with sufficient node count or clone count.
2.  **Construct global network from the filtered data and perform
    clustering.** Combine the filtered data from all samples into a
    single global network. Perform cluster analysis and assign
    membership to the global clusters. Assess sample representation
    within the clusters to identify public clusters.

## Simulate Data for Demonstration

We simulate some toy data for demonstration.

We simulate a total of 30 samples, each containing 30 observations.

Some sequences are simulated with a tendency to appear in relatively few
samples, while others are simulated with a tendency to appear in many
samples.

Each sample is saved in a separate file using the .rds file format. The
files are named “`Sample1.rds`”, “`Sample2.rds`”, etc. The file path of
their directory is saved to the R environment variable
`dir_input_samples` for later reference.

``` r
# Use temp dir
data_dir <- tempdir()

# Directory to store input files
dir_input_samples <- file.path(data_dir, "input_samples")
dir.create(dir_input_samples, showWarnings = FALSE)

samples <- 30
sample_size <- 30 # (seqs per sample)          

base_seqs <- c(
  "CASSIEGQLSTDTQYF", "CASSEEGQLSTDTQYF", "CASSSVETQYF",
  "CASSPEGQLSTDTQYF", "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF",
  "CASSLTSGYNEQFF", "CASSETGYNEQFF", "CASSLTGGNEQFF", "CASSYLTGYNEQFF",
  "CASSLTGNEQFF", "CASSLNGYNEQFF", "CASSFPWDGYGYTF", "CASTLARQGGELFF",
  "CASTLSRQGGELFF", "CSVELLPTGPLETSYNEQFF", "CSVELLPTGPSETSYNEQFF",
  "CVELLPTGPSETSYNEQFF", "CASLAGGRTQETQYF", "CASRLAGGRTQETQYF",
  "CASSLAGGRTETQYF", "CASSLAGGRTQETQYF", "CASSRLAGGRTQETQYF",
  "CASQYGGGNQPQHF", "CASSLGGGNQPQHF", "CASSNGGGNQPQHF", "CASSYGGGGNQPQHF",
  "CASSYGGGQPQHF", "CASSYKGGNQPQHF", "CASSYTGGGNQPQHF", 
  "CAWSSQETQYF", "CASSSPETQYF", "CASSGAYEQYF", "CSVDLGKGNNEQFF") 

# relative generation probabilities
pgen <- cbind(
  stats::toeplitz(0.6^(0:(sample_size - 1))),
  matrix(1, nrow = samples, ncol = length(base_seqs) - samples))

# Simulate the data
library(NAIR)
foo <- simulateToyData(    
  samples = samples, sample_size = sample_size,
  prefix_length = 1, prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, samples), rep(0, samples)),
  affixes = base_seqs, affix_probs = pgen, num_edits = 0,
  output_dir = dir_input_samples, no_return = FALSE)
```

The first few rows of the data for the first sample appear as follows:

``` r
# View first few rows of data for sample 1
head(readRDS(file.path(dir_input_samples, "Sample1.rds")))
#>           CloneSeq CloneFrequency CloneCount SampleID
#> 1 CASSIEGQLSTDTQYF     0.02606559       2832  Sample1
#> 2 CASSEEGQLSTDTQYF     0.03718396       4040  Sample1
#> 3      CASSSPETQYF     0.03182726       3458  Sample1
#> 4 CASSIEGQLSTDTQYF     0.04615781       5015  Sample1
#> 5      CAWSSQETQYF     0.06006498       6526  Sample1
#> 6 CASSEEGQLSTDTQYF     0.03363123       3654  Sample1
```

## 1. Filter Clusters Within Each Sample

First, we use the `findPublicClusters` function to perform network
analysis on each sample individually and select clusters based on node
count and clone count.

Before explaining its usage in detail, we briefly describe what the
function does. For each sample, the repertoire network is constructed
and cluster analysis is used to partition the network into clusters. The
clusters are then filtered according to node count and clone count based
on user-specified criteria. The node-level and cluster-level meta data
for the clusters that remain after filtering are saved as files to be
used as inputs for [step 2](#2.-global-network-of-public-clusters).

### 1.1 Filter Settings

The `findPublicClusters` function has several parameters that control
the criteria used to filter the nodes and clusters in each sample. These
arguments are presented below.

#### 1.1.1 Top $n$ Clusters

Within each sample, the clusters are ranked by node count. The top $n$
highest-ranking clusters (those with the greatest node count) within
each sample are automatically retained. The default value of $n$ is 20.
A different value of $n$ can be specified using the `top_n_clusters`
argument.

If more than one cluster is tied for $n$-th place in the ranking by node
count, only one such cluster will be included in the top $n$ clusters.
This ensures that no more than $n$ clusters are selected from each
sample based on this criterion.

If fewer than $n$ clusters are present in the network for a sample, then
all of the clusters will be retained.

#### 1.1.2 Minimum Node Count

In addition to retaining the top $n$ clusters from each sample, clusters
that contain a sufficient number of nodes will also be retained. By
default, any cluster containing at least ten nodes will be retained.
This value can be adjusted using the `min_node_count` argument. For
example, setting `min_clone_count = 30` will retain all clusters
containing at least 30 nodes.

#### 1.1.3 Minimum Clone Count

In addition to the clusters retained based on node count, clusters with
a sufficient aggregate clone count will also be retained. The aggregate
clone count of a cluster is the sum of the clone counts across all nodes
(clones) in the cluster.

By default, any cluster with an aggregate clone count of at least 100
will be retained. This value can be adjusted using the `min_clone_count`
argument. For example, setting `min_clone_count = 500` will retain all
clusters with an aggregate clone count of at least 500.

#### 1.1.4 Sequence Length

When building the network for each sample, only clones whose receptor
sequences are at least three characters in length will be included in
the network. This minimum value for sequence length can be adjusted by
setting the `min_seq_length` argument to a different value. Setting the
value to `NULL` bypasses this check.

#### 1.1.5 Sequence Content

When building the network for each sample, all clones whose receptor
sequences contain characters `*` or `_` will be omitted from the
network. This can be changed using the `drop_matches` argument, which
accommodates a character string or regular expression specifying the
pattern of content to search for. The content of each clone’s sequence
is checked for a match to this pattern using the `grep` function from
base R. If a match is found, the clone is omitted from the network.
Setting the value to `NULL` bypasses this check.

For details on how the pattern matching is performed, please refer to
the base R documentation files for `regex` and `grep`.

### 1.2 Input File List

The main argument of the `findPublicClusters` function is the
`file_list` argument, which accepts a vector containing file paths. Each
path corresponds to a distinct AIRR-Seq data file representing an
individual sample.

Below, we prepare the vector `input_files` to be provided to the
`file_list` argument of `findPublicClusters`:

``` r
# input files for step 1 (one per sample)
input_files <- file.path(dir_input_samples, paste0("Sample", 1:samples, ".rds"))
head(input_files)
#> [1] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample1.rds"
#> [2] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample2.rds"
#> [3] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample3.rds"
#> [4] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample4.rds"
#> [5] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample5.rds"
#> [6] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/input_samples/Sample6.rds"
```

### 1.3 Input Type

The file format of the input files for `findPublicClusters` is specified
using the `input_type` parameter. The supported formats include `"rds"`,
`"rda"`, `"csv"`, as well as files that can be read using the
`read.table` function, such as `"tsv"` and `"txt"`.

For text formats such as `"csv"`, `"tsv"` and `"txt"`, users can specify
the separation option by utilizing the `sep` argument. The default
setting `sep = ""` accommodates all forms of white space, i.e., one or
more spaces, tabs, newlines or carriage returns. In addition, it is
important to note that the first line of the data is assumed to be the
header by default. To disable this behavior and treat the first line as
data, users must set the `header` parameter to `FALSE`.

Our samples are stored in .rds files, so we use `input_type = "rds"`.

### 1.4 Specifying the Sequence Column

The `seq_col` argument is used to specify the column containing the
clone sequences in the input data for each sample. The argument accepts
either the column name or column index.

In our simulated data, the column containing the clone sequences is
named `CloneSeq`.

### 1.5 Assigning Custom Sample IDs

The sample ID for each clone will be included in the filtered data saved
by `findPublicClusters`. This will allow us to distinguish the sample of
origin for each clone after we combine the filtered data from all
samples in [step 2](#2.-global-network-of-public-clusters).

By default, the samples are labeled numerically according to the order
they appear in `file_list`. The `sample_ids` argument allows for custom
sample IDs to be assigned, if desired. The argument accepts a vector of
the same length as `file_list`. Each entry of `sample_ids` is assigned
as a sample ID to the sample in the corresponding entry of `file_list`.

### 1.6 Customizing the Network Analysis

Most of the arguments from
[`buildRepSeqNetwork`](buildRepSeqNetwork.html) used to customize the
network analysis can also be used in `findPublicClusters` to customize
the network analysis for each sample. We mention some of the important
ones below.

Note that the arguments `node_stats`, `stats_to_include` and
`cluster_stats` are not available, since `findPublicClusters` always
computes a fixed set of node-level and cluster-level properties.

#### 1.6.1 Distance Function and Distance Cutoff

By default, two nodes within a sample’s network are joined by an edge if
the Hamming distance between their sequences is at most 1.

The Levenshtein distance can be used instead of the Hamming distance by
setting `dist_type = "lev"`.

The maximum distance for two nodes to be joined by an edge can be
changed using the `dist_cutoff` argument.

#### 1.6.2 Clustering Algorithm

By default, clustering is performed using the `cluster_fast_greedy`
algorithm from the `igraph` package. Other clustering algorithms from
the `igraph` package can be used instead of the default algorithm. The
algorithm is specified using the `cluster_fun` argument, which accepts
one of the following functions:

- `cluster_fast_greedy`
- `cluster_edge_betweenness`
- `cluster_fluid_communities`
- `cluster_infomap`
- `cluster_label_prop`
- `cluster_leading_eigen`
- `cluster_leiden`
- `cluster_louvain`
- `cluster_optimal`
- `cluster_spinglass`
- `cluster_walktrap`

For example, setting `cluster_fun = cluster_leiden` performs clustering
using the `cluster_leiden` algorithm.

For more information about a particular algorithm, users can refer to
its help documentation file. For example, the command
`?cluster_fast_greedy` loads the documentation file for the
`cluster_fast_greedy` (the default) algorithm, assuming the `NAIR`
package has been loaded (e.g., using `library(NAIR)`).

### 1.7 Output Settings

The `findPublicClusters` function does not return any direct output.
Instead, it saves the network meta data for the selected clusters to
files that will be used as inputs in step 2.

For each sample, `findPublicClusters` saves two files, one containing
the node-specific network meta data and the other containing the
cluster-specific network meta data. The data for each sample is saved
after filtering the clusters, and thus includes only the clusters
selected based on the filtering criteria detailed in section
[1.1](#1.1-filter-settings).

The file path for the output directory is specified using the
`output_dir` argument. The output directory specified by the
`output_dir` argument will be created if it does not already exist.

Within the output directory, two subdirectories are created. One
subdirectory, named `node_meta_data`, contains the node-level data files
for each sample, and the other, named `cluster_meta_data`, contains the
cluster-level data files. Within each of these two subdirectories, the
file for each sample is named according to its sample ID as specified by
the [`sample_ids` argument](#1.5-assigning-custom-sample-ids).

By default, each file is saved as a RDS file. This can be changed using
the `output_type` argument. Other valid options include `"rda"` and
`"csv"`.

#### 1.7.1 Saving Unfiltered Network Data

By default, the `findPublicClusters` function saves the network meta
data only for the clusters selected based on the filtering criteria
detailed in section [1.1](#1.1-filter-settings).

If desired, the network meta data for each sample’s full network can
also be saved prior to filtering the clusters. This is done by by
providing a file path to the `output_dir_unfiltered` argument, which
specifies a separate output directory for the full network data. Each
sample’s full network data prior to filtering the clusters will then be
saved to the directory specified by the `output_dir_unfiltered`
argument. This data is saved separately from, and in addition to, the
[default data that is saved after filtering the
clusters](#1.7-output-settings).

Note that the sequence-based filter settings specified by
[`min_seq_length`](#1.1.4-sequence-length) and
[`drop_chars`](#1.1.5-sequence-content) still apply to each sample’s
full network, since the network is only constructed after applying these
filters. The full, pre-filtered network refers to the network that
contains all of the sample’s clusters, i.e., before the public clusters
are identified and other clusters are removed.

The full network data for each sample includes node-level meta data,
cluster-level meta data, the `igraph` edge list, the network adjacency
matrix, as well as any plots generated ([see next
subsection](#1.8-(optional)-network-visualization-per-sample)). By
default, the R objects for these files are saved into a single RDS file
whose file name (excluding the .rds file extension) is the [sample
ID](#1.5-assigning-custom-sample-ids). The RDA file format can be used
instead of RDS by setting `output_type_unfiltered = "rda"`. If desired,
each R object can be saved to a separate file by setting
`output_type_unfiltered = "individual"`. This saves each object
according to the default file format used by
[`buildRepSeqNetwork`](buildRepSeqNetwork.html). When saving objects
individually, the [sample ID](#1.5-assigning-custom-sample-ids) is used
as the common file name prefix for the files from each sample.

### 1.8 (Optional) Network Visualization Per Sample

By default, `findPublicClusters` does not produce visual plots when
constructing the network for each sample. Instead, visualization occurs
after combining the data from all samples into a single network in [step
2](#2.-global-network-of-public-clusters).

If desired, visual plots of the full network for each sample ([prior to
filtering the clusters](#1.7.1-saving-unfiltered-network-data)) can be
produced by setting `plots = TRUE`. The plots can be printed in R by
additionally setting `print_plots = TRUE`. Furthermore, if the user has
specified to [save the full network data for each sample (prior to
filtering the clusters)](#1.7.1-saving-unfiltered-network-data) by
providing a file path to the `output_type_unfiltered` argument, these
plots will be saved along with the rest of the full network data for
each sample.

Note that the plots are not saved along with the [default data that is
saved after filtering the clusters](#1.7-output-settings). Therefore, if
no file path is provided to the `output_type_unfiltered` argument, the
plots will not be saved at all. If, in addition, the argument
`print_plots` is set to `FALSE` (the default), setting `plots = TRUE`
will have no effect.

If the user sets `plots` to `TRUE`, then by default, the plot for each
sample will color the nodes according to their cluster membership. If
desired, a different variable can be used to color the nodes. This is
done by specifying the column for the variable to the `color_nodes_by`
argument, which accepts a column name or column index.

Other arguments accepted by `buildRepSeqNetwork` to customize the
visualization can also be used when calling `findPublicClusters`.

### 1.9 Execution and Output

We execute the `findPublicClusters` function using the input we prepared
earlier for the `file_list` argument:

``` r
# 1. Filter Clusters Within Each Sample
dir_filtered_samples <- file.path(data_dir, "filtered_samples")
findPublicClusters(
  file_list = input_files, input_type = "rds",
  sample_ids = paste0("Sample", 1:samples),
  seq_col = "CloneSeq", count_col = "CloneCount",
  min_seq_length = NULL, drop_matches = NULL,
  top_n_clusters = 3, min_node_count = 5, min_clone_count = 15000,
  output_dir = dir_filtered_samples)
#> <<< Beginning search for public clusters >>>
#> Processing sample 1 of 30: Sample1
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 30 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (18 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 2 of 30: Sample2
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 6 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (22 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 3 of 30: Sample3
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 29 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 6 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (20 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 4 of 30: Sample4
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 6 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (18 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 5 of 30: Sample5
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 27 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 6 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (22 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 6 of 30: Sample6
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 29 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 9 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (15 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 7 of 30: Sample7
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (13 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 8 of 30: Sample8
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 9 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (12 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 9 of 30: Sample9
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (13 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 10 of 30: Sample10
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 9 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (15 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 11 of 30: Sample11
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (18 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 12 of 30: Sample12
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 29 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 10 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (12 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 13 of 30: Sample13
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 27 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (18 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 14 of 30: Sample14
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (17 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 15 of 30: Sample15
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 29 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (14 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 16 of 30: Sample16
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 6 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (21 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 17 of 30: Sample17
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 29 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (22 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 18 of 30: Sample18
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 9 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (14 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 19 of 30: Sample19
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (17 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 20 of 30: Sample20
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (15 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 21 of 30: Sample21
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (19 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 22 of 30: Sample22
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (16 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 23 of 30: Sample23
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 27 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (19 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 24 of 30: Sample24
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 27 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 8 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (18 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 25 of 30: Sample25
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 9 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (17 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 26 of 30: Sample26
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (17 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 27 of 30: Sample27
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 25 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 5 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (18 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 28 of 30: Sample28
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 26 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (16 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 29 of 30: Sample29
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 29 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 7 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 4 clusters (22 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> Processing sample 30 of 30: Sample30
#> Input data contains 30 rows.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 28 nodes (after removing isolated nodes).
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 6 clusters in the network... Done.
#> >>> Filtering clusters in the current sample... Done.
#> * 3 clusters (19 nodes) remain. Saving results... Done.
#> ----------------------------------------------------------------------
#> All samples complete. Filtered data is located in the following directory:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/filtered_samples
```

## 2. Global Network of Public Clusters

Next, we use the `buildPublicClusterNetwork` function to combine the
filtered data from all samples into a single global network and perform
clustering analysis.

### 2.1 Input File List

The files created by `findPublicClusters` in the previous step contain
the filtered data for each sample. As [detailed
earlier](#1.7-output-settings), these files are located in two separate
subdirectories, one containing the files for the node-level meta data
and the other containing the files for the cluster-level meta data.

For this step, we require only the files containing the node-level meta
data. These files are provided to `buildPublicClusterNetwork` by
supplying a character vector of file paths to the `file_list` argument.
We create this vector below.

``` r
# Node-level meta data for each sample's filtered clusters
dir_filtered_samples_node <- file.path(dir_filtered_samples, "node_meta_data")
files_filtered_samples_node <- list.files(dir_filtered_samples_node, 
                                          full.names = TRUE)
head(files_filtered_samples_node)
#> [1] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/filtered_samples/node_meta_data/Sample1.rds" 
#> [2] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/filtered_samples/node_meta_data/Sample10.rds"
#> [3] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/filtered_samples/node_meta_data/Sample11.rds"
#> [4] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/filtered_samples/node_meta_data/Sample12.rds"
#> [5] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/filtered_samples/node_meta_data/Sample13.rds"
#> [6] "C:\\Users\\Brian\\AppData\\Local\\Temp\\RtmpsFKIdb/filtered_samples/node_meta_data/Sample14.rds"
```

### 2.2 Customization of Network Analysis

`buildPublicClusterNetwork` uses the same arguments as
[`buildRepSeqNetwork`](buildRepSeqNetwork.html) for customizing the
network analysis. We mention some of the important ones below.

#### 2.2.1 Distance Function and Distance Cutoff

By default, two nodes within the global network are joined by an edge if
the Hamming distance between their sequences is at most 1.

The Levenshtein distance can be used instead of the Hamming distance by
setting `dist_type = "lev"`.

The maximum distance for two nodes to be joined by an edge can be
changed using the `dist_cutoff` argument.

#### 2.2.2 Clustering Algorithm

After constructing the global network, `buildPublicClusterNetwork`
performs cluster analysis on the network nodes, partitioning the global
network graph into densely-connected subgraphs. These global clusters
can contain nodes from different samples.

By default, the clustering is performed using the `cluster_fast_greedy`
algorithm. The clustering algorithm can be changed using the
`cluster_fun` argument [as described
earlier](#1.6.2-clustering-algorithm).

### 2.3 Customization of Visual Plot

By default, the network graph plot produced by
`buildPublicClusterNetwork` colors the nodes according to sample ID.
This can assist the user in identifying the public clusters. If desired,
a different variable can be used to color the nodes. This is done by
specifying the column for the variable to the `color_nodes_by` argument,
which accepts a column name or column index.

Other arguments accepted by `buildRepSeqNetwork` to customize the
visualization can also be used when calling `buildPublicClusterNetwork`.

### 2.4 Output Settings

The output returned by `buildPublicClusterNetwork` follows the same
format as the output of [`buildRepSeqNetwork`](buildRepSeqNetwork.html).
The function returns a list containing the node-level and cluster-level
meta data for the global network, as well as any plots generated, in
addition to the network adjacency matrix and the `igraph` network edge
list.

By default, a subdirectory named `"public_clusters"` is created within
the current working directory, and the contents of the list returned by
`buildPublicClusterNetwork` are saved to this subdirectory. Each list
element is saved as an individual file. The file formats are the same
default file formats used by
[`buildRepSeqNetwork`](buildRepSeqNetwork.html). In particular, the
node-level and cluster-level meta data are saved as csv files.

Alternatively, the user can save the entire output list to a single
compressed rds or rda file by setting `output_type = "rds"` or
`output_type = "rda"`, respectively.

By default, all files saved share the common file name prefix
`"PublicClusterNetwork"`. This common file name prefix can be set to a
different value by supplying a character string to the `output_name`
argument.

The output can be saved to a different directory by providing a file
path to the `output_dir` argument.

The user can also specify `output_dir = NULL` in order to prevent the
output from being saved.

### 2.5 Execution and Output

We execute the `buildPublicClusterNetwork` function using the input we
prepared earlier for the `file_list` argument:

``` r
dir_out <- file.path(data_dir, "public_clusters")

# Collect clones from all public clusters and perform network analysis
public_clusters <- buildPublicClusterNetwork(
  file_list = files_filtered_samples_node,
  seq_col = "CloneSeq", count_col = "CloneCount",
  size_nodes_by = 1,
  output_dir = dir_out)
#> Building network of public clusters:
#> Input data contains 517 rows.
#> Removing sequences with length fewer than 3 characters... Done. 517 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 517 nodes.
#> Computing cluster membership within the network... Done.
#> Computing node-level network statistics... Done.
#> Computing statistics for the 20 clusters in the network... Done.
#> Generating graph plot with nodes colored by SampleID...
#>  Done.
#> Node-level meta-data saved to file:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/public_clusters/PublicClusterNetwork_NodeMetadata.csv
#> Cluster-level meta-data saved to file:
#>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/public_clusters/PublicClusterNetwork_ClusterMetadata.csv
```

<img src="man/figures/README-unnamed-chunk-49-1.png" width="100%" style="display: block; margin: auto;" />

    #> Network graph plots saved to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/public_clusters/PublicClusterNetwork.pdf
    #> Network igraph saved in edgelist format to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/public_clusters/PublicClusterNetwork_EdgeList.txt
    #> Adjacency matrix saved to file:
    #>   C:\Users\Brian\AppData\Local\Temp\RtmpsFKIdb/public_clusters/PublicClusterNetwork_AdjacencyMatrix.mtx

The elements of the output list returned by `buildPublicClusterNetwork`
have the following names:

``` r
# output returned by buildPublicClusterNetwork
names(public_clusters)
#> [1] "igraph"           "adjacency_matrix" "node_data"        "cluster_data"    
#> [5] "plots"
```

The default files saved by `buildPublicClusterNetwork` have the
following names:

``` r
# Files saved by buildPublicClusterNetwork
list.files(dir_out)
#> [1] "PublicClusterNetwork.pdf"                
#> [2] "PublicClusterNetwork_AdjacencyMatrix.mtx"
#> [3] "PublicClusterNetwork_ClusterMetadata.csv"
#> [4] "PublicClusterNetwork_EdgeList.txt"       
#> [5] "PublicClusterNetwork_NodeMetadata.csv"
```

#### 2.5.1 Node-Level Meta Data

The `node_data` data frame contained in the output list contains the
following variables:

``` r
# variables in the node-level meta data
names(public_clusters$node_data)
#>  [1] "CloneSeq"                           "CloneFrequency"                    
#>  [3] "CloneCount"                         "SampleID"                          
#>  [5] "SampleLevelNetworkDegree"           "ClusterIDInSample"                 
#>  [7] "SampleLevelTransitivity"            "PublicCloseness"                   
#>  [9] "SampleLevelCentralityByCloseness"   "SampleLevelEigenCentrality"        
#> [11] "SampleLevelCentralityByEigen"       "SampleLevelBetweenness"            
#> [13] "SampleLevelCentralityByBetweenness" "SampleLevelAuthorityScore"         
#> [15] "SampleLevelCoreness"                "SampleLevelPageRank"               
#> [17] "PublicNetworkDegree"                "ClusterIDPublic"                   
#> [19] "PublicTransitivity"                 "PublicCentralityByCloseness"       
#> [21] "PublicEigenCentrality"              "PublicCentralityByEigen"           
#> [23] "PublicBetweenness"                  "PublicCentralityByBetweenness"     
#> [25] "PublicAuthorityScore"               "PublicCoreness"                    
#> [27] "PublicPageRank"
```

Notice that by default, all variables that were present in each sample’s
original data, such as `"CloneFrequency"` and `"CloneCount"`, are
automatically carried over into this data.

Many node-level network properties are present in the data. Some pertain
to the network for the individual sample from which each node
originated, i.e., the networks constructed in [step
1](#1.-filter-clusters-within-each-sample). This includes all variables
that begin with `"SampleLevel"`, such as `SampleLevelNetworkDegree`, as
well as the variable `ClusterIDInSample`.

The variables that begin with `"Public"`, such as `PublicNetworkDegree`,
pertain to the global network. In particular, `ClusterIDPublic`
indicates the ID of the global cluster to which each node belongs.

#### 2.5.2 Cluster-Level Meta Data

The `cluster_data` data frame contained in the output list contains the
following variables:

``` r
# variables in the node-level meta data
names(public_clusters$cluster_data)
#>  [1] "cluster_id"                  "node_count"                 
#>  [3] "mean_seq_length"             "mean_degree"                
#>  [5] "max_degree"                  "seq_w_max_degree"           
#>  [7] "agg_count"                   "max_count"                  
#>  [9] "seq_w_max_count"             "diameter_length"            
#> [11] "global_transitivity"         "assortativity"              
#> [13] "edge_density"                "degree_centrality_index"    
#> [15] "closeness_centrality_index"  "eigen_centrality_index"     
#> [17] "eigen_centrality_eigenvalue"
```

Each row corresponds to a cluster in the global network, and each
variable corresponds to a cluster-level property.

### 2.6 Labeling the Global Clusters

In the plot, there appear to be six clusters with more than 3 samples
represented.

In order to reference these clusters within the data, we can [label the
six largest clusters in the plot with their cluster
IDs](network_visualization.html#labeling-clusters) using the
`addClusterLabels` function. Note that within the node-level metadata,
the global cluster ID is stored in the variable `ClusterIDPublic`, so we
must provide this column name to the `cluster_id_col` argument.

``` r
# Modify plot to add labels to the clusters
public_clusters$plots[[1]] <- 
  addClusterLabels(
    plot = public_clusters$plots[[1]],
    net = public_clusters,
    top_n_clusters = 6,
    cluster_id_col = "ClusterIDPublic",
    size = 7
  )

# View modified plot
public_clusters$plots[[1]]
#> Warning: Removed 511 rows containing missing values (`geom_text()`).
```

<img src="man/figures/README-unnamed-chunk-54-1.png" width="100%" style="display: block; margin: auto;" />

### 2.7 Focusing on Individual Clusters of Interest

If we wish to focus on a particular cluster of interest within the
global network, we can build a network exclusively using the clones from
that cluster. This is accomplished using the `buildRepSeqNetwork`
function, where we subset our data according to the value of the
`ClusterIDPublic` column, which contains the global cluster IDs.

Below, we focus on the first cluster, which in this case is also the
largest cluster by node count. In the plot, we color each node according
to the receptor sequence of its corresponding clone.

``` r
# focus on cluster 1
buildRepSeqNetwork(
  data = 
    public_clusters$node_data[public_clusters$node_data$ClusterIDPublic == 1, ],
  seq_col = "CloneSeq", 
  color_nodes_by = "CloneSeq", color_scheme = "plasma", 
  size_nodes_by = 3, 
  output_dir = NULL, output_name = "Cluster 1")
#> Input data contains 96 rows.
#> Removing sequences with length fewer than 3 characters... Done. 96 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 96 nodes (after removing isolated nodes).
#> Generating graph plot with nodes colored by CloneSeq...
```

<img src="man/figures/README-unnamed-chunk-55-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.

Below, we do the same for the sixth largest cluster:

``` r
# focus on cluster 6
buildRepSeqNetwork(
  data = 
    public_clusters$node_data[public_clusters$node_data$ClusterIDPublic == 6, ],
  seq_col = "CloneSeq", 
  color_nodes_by = "CloneSeq", color_scheme = "plasma", 
  size_nodes_by = 3, 
  output_dir = NULL, output_name = "Cluster 6")
#> Input data contains 27 rows.
#> Removing sequences with length fewer than 3 characters... Done. 27 rows remaining.
#> Computing network edges based on a max hamming distance of 1... Done.
#> Network contains 27 nodes (after removing isolated nodes).
#> Generating graph plot with nodes colored by CloneSeq...
```

<img src="man/figures/README-unnamed-chunk-56-1.png" width="100%" style="display: block; margin: auto;" />

    #>  Done.
