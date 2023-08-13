

# Changes in NAIR version 0.0.9044 (8/13/2023)

## Package 

* Tests
    * Tests added for `hamDistBounded()`, `levDistBounded()`, `sparseAdjacencyMatFromSeqs()`.
    * Tests added for low-level argument checks.
* DESCRIPTION
    * Added R minimum version requirement 3.1.0 to `Depends`, since version 3.0.2 or greater is needed to require specific minimum versions of RcppArmadillo and Rcpp in the `LinkingTo` field (requiring 3.1.0 since CRAN advises against requiring R versions that don't have 0 as the third value).
* Compile Flags in Makevars
    * Re-enabled OpenMP support for Unix-based systems.

## Documentation and Vignettes

* URLs in documentation files and vignettes have been curated to conform with CRAN's policies. Specifically, any URL that redirected to another URL has been replaced by the target of the redirection. Links to package CRAN pages are now in canonical form.
* All function reference files now run their examples when the package is built or checked. Some examples have been expanded.
* Examples and vignette code now remove files and directories created in the temporary directory for the current R session.
* Added a documentation file for the package itself (`NAIR-package`).
* All function reference files now have a Value section, including functions that do not return a value.

## Functions

* Argument checks have been expanded to encompass most user-facing package functions.
* `addClusterMembership()`
    * The `fun` argument is now passed to `match.fun()` before being called. This allows users to specify clustering algorithms using the syntax, e.g., `cluster_fun = "cluster_walktrap"` in addition to the previously-accepted `cluster_fun = cluster_walktrap`.
* `findAssociatedClones()`
    * Upon completion, removes files and directories created in the temporary directory for the current R session.
* Internal package C++ functions `hamAdjacencyMatSparse()` and `levAdjacencyMatSparse()` have an additional argument specifying the directory in which to write the file `col_ids.txt` (previously it was written to the current working directory of the R session at the time they were called).
* `sparseAdjacencyMatFromSeqs()`
    * Now uses the temporary directory for the current R session when calling internal package C++ functions `hamAdjacencyMatSparse()` and `levAdjacencyMatSparse()`, rather than the current working directory.



# Changes in NAIR version 0.0.9043 (8/7/2023)

## Vignettes

* Function names in vignettes have been reformatted so that when pkgdown builds the package webpage articles from the vignettes, function names will link to the webpages with their documentation files.

## Readme

* Package README file updated

## Website

* Custom index added for the reference topics. Topics are now organized into named sections, dramatically improving the ability to find particular topics or functions of interest.

## Package

* `packageStartupMessage()` added to `.onAttach()`: When loaded, the package will provide a welcome message with instructions for getting started.


# Changes in NAIR version 0.0.9042 (8/6/2023)

## Package Functions

* `findAssociatedClones`
    * The variable `SampleID` created in the output data is now forced to be of type character. Previously, values from the argument `sample_ids` were sometimes unintentionally converted from character to numeric, such as the default values in `sample_ids`, which were `"1"`, `"2"`, etc. This was causing these variables to be treated as continuous variables when used to color nodes in the network graph plot, which resulted in their color scales being depicted in the wrong format in the plot legend.
    * The `sample_ids` argument is now coerced to a character vector. This prevents an error when saving the output that occurred when `sample_ids` used numeric values.
    * Default value of `sample_ids` now has entries `"Sample1"`, `"Sample2"`, etc., instead of `"1"`, `"2"`, etc. 
* `findPublicClusters`
    * The variables `SampleID`, `SubjectID` and `GroupID` created in the output data are now forced to be of type character. Previously, values from the arguments `sample_ids`, `subject_ids` and `group_ids` were sometimes unintentionally converted from character to numeric, such as the default values in `sample_ids`, which were `"1"`, `"2"`, etc. This was causing these variables to be treated as continuous variables when used to color nodes in the network graph plot, which resulted in their color scales being depicted in the wrong format in the plot legend.
    * The `sample_ids` argument is now coerced to a character vector. This prevents an error when saving the output that occurred when `sample_ids` used numeric values (which was the previous default!).
    * Default value of `sample_ids` now has entries `"Sample1"`, `"Sample2"`, etc., instead of `1`, `2`, etc. 
* `buildPublicClusterNetwork`
    * Argument `plot_title` added with default value `"Global Network of Public Clusters"`. Previously this argument was passed to `buildRepSeqNetwork` through the ellipses `...` argument, and thus used a default value of `"auto"`, which resulted in the default plot title being the value of the `output_name` argument, which is `"PublicClusterNetwork"` by default. 
* `plotNetworkGraph`
    * Added arguments `pdf_width` and `pdf_height` for adjusting the dimensions of the pdf when saving this function's output directly using the `outfile` argument. Other package functions use `saveNetworkPlots` for saving plots created using `plotNetworkGraph`, so the absence of these arguments in the `plotNetworkGraph` function had gone unnoticed previously. But since the function has an option to save the output directly to pdf using the `outfile` argument, it is only appropriate to also provide control over the pdf dimensions.
* `kmeansAtchley`
    * Default values for `amino_col` and `sample_col` arguments removed as the previous defaults are no longer useful. They were originally designed based on a previous version of the associated clusters workflow.
    * Default filenames for the pdfs of the heatmaps have been changed to `"atchley_kmeans_TCR_fraction_per_cluster.pdf"` and `"atchley_kmeans_correlation_heatmap.pdf"`. The previous values `"atchley_kmeans_cluster_relative_size_profiles_by_sample.pdf"` and `"atchley_kmeans_corr_in_cluster_size_profile_between_samples.pdf"` were longer and potentially more confusing in their meaning.
    
## Vignettes

* `Searching for Associated TCR/BCR Clusters`
    * Completion of major revisions
* `Searching for Public Clusters`
    * Completion of major revisions
* `buildRepSeqNetwork`
    * Minor content revisions
* `Network Visualization`
    * Restructured sections and minor content revisions
    
## Documentation

* All package documentation files have been revised and updated.


# Changes in NAIR version 0.0.9041 (8/2/2023)

## Package Functions

* `buildAssociatedClusterNetwork`
    * Now performs clustering analysis to obtain cluster membership ID even if the user manually specifies not to compute cluster stats or the `cluster_id` network property. This is because performing clustering and obtaining the cluster membership is a primary purpose of this function. It is still desirable for the user to be able to prevent other node-level properties as well as cluster-level properties from being computed if desired, but now doing so will not interfere with the function accomplishing its purpose. 
* `findPublicClusters`
    * Fixed a bug that was causing the sample-level network node property `SampleLevelCloseness` to be left named as `closeness` in the data frames for the filtered node-level data. This bug was in turn causing this property to be overwritten by the global network node property `PublicCloseness` when calling `buildPublicClusterNetwork`. 
    
## Vignettes

* `Searching for Associated TCR/BCR Clusters`
    * Completion of major revisions
* `Searching for Public Clusters`
    * Completion of major revisions
* `buildRepSeqNetwork`
    * Minor content revisions
* `Network Visualization`
    * Restructured sections and minor content revisions

# Changes in NAIR version 0.0.9040

## Package Functions

* `buildAssociatedClusterNetwork`
    * Default value of `data_symbols` argument changed from `NULL` to `"data"` in order to match the output format of `findAssociatedClones` when `findAssociatedClones` is called with `output_type = "rda"`. Note this change only affects the case when `buildAssociatedClusterNetwork` is called with `input_type = "rda"`
    * Now performs clustering analysis to obtain cluster membership ID even if the user manually specifies not to compute cluster stats or the `cluster_id` network property. This is because performing clustering and obtaining the cluster membership is a primary purpose of this function. It is still desirable for the user to be able to prevent other node-level properties as well as cluster-level properties from being computed if desired, but now doing so will not interfere with the function accomplishing its purpose. 
    
## Vignettes

* `buildRepSeqNetwork`
    * For both `buildRepSeqNetwork` and `saveNetwork`, the vignette now specifies the R environment variable name for the output list when it is saved to an Rdata file using `output_type = "rda"`.
* `Searching for Associated TCR/BCR Clusters`
    * Content added to more completely explain certain behavior and arguments that were not previously covered.
    
## Package Webpage

* GitHub Actions workflow added to automate publication of future webpage updates 
* URL for GitHub Pages hosted website added to URL field of DESCRIPTION file and to pkgdown.yaml

# Changes in NAIR version 0.0.9039 (7/29/2023)

## Package Metadata

* Updated authors in DESCRIPTION
* Updated journal article citation in Readme and documentation files

## Vignettes

* The term "AIRR-Seq" is now spelled out in full as "Adaptive Immune Receptor Repertoire Sequencing" prior to the first instance of its abbreviation in each vignette in which it appears
* `buildRepSeqNetwork`
    * Meaning of cluster membership ID and corresponding `cluster_id` network property are now more clearly explained
* `Searching for Associated TCR/BCR Clusters`
    * Header levels and document structure updated
    * Further revisions/additions, including fixing the remaining broken links, are forthcoming
    
## Documentation

* References/Authors updated across all documentation files
* Many documentation files have been revised and updated to use wording that is more clear, accurate and consistent with the language used in package vignettes. Updates to the remaining documentation files are forthcoming.  


# Changes in NAIR version 0.0.9038 (7/20/2023)

## Bug Fixes

* Fixed a bug in `filterInputData` that raised an error when the `count_col` and `subset_cols` arguments were both non-null

## Functions

* When calling `plotNetworkGraph` directly with a vector provided to `color_nodes_by` and with `color_title = "auto"` (the default), the function will attempt to use the name of the vector for the color legend title. A similar change applies with respect to the arguments `size_nodes_by` and `size_title`.
* `buildPublicClusterNetwork` arguments `node_stats`, `stats_to_include` and `cluster_stats` are now deprecated and do nothing. All node-level and cluster-level network properties are now automatically computed. The arguments remain in order to maintain backwards compatibility with user code, but raise a warning notifying the user of their deprecated state when a non-null value is provided.
* Functions for clustering algorithms imported from the `igraph` package (such as `cluster_fast_greedy`) are now exported in the package NAMESPACE file so that they are available to users. These functions can now be used as inputs to the `cluster_fun` argument of various `NAIR` package functions without the need to use the `igraph::` prefix.
* Documentation file added for the above re-exported functions
    
## Vignettes and Documentation

* `Utility Functions` vignette (formerly titled `Downstream Analysis`) removed. Its content has been absorbed into the `buildRepSeqNetwork` and `Network Visualization` vignettes
* Revisions and content additions to the following vignettes:
    * `buildRepSeqNetwork`
    * `Network Visualization`
    * `Searching for Public Clusters`
* The help file for `plotNetworkGraph` now recommends that users prefer the higher-level function `generateNetworkGraphPlots` over `plotNetworkGraph`, since the former has arguments that behave identically to those of `buildRepSeqNetwork` and supports generation of multiple plots. `plotNetworkGraph` is called by `generateNetworkGraphPlots`, so users should have no need to call `plotNetworkGraph` directly. However, `plotNetworkGraph` remains as an exported function available to the user in order to maintain backwards compatibility with user code. 


# Changes in NAIR version 0.0.9037 (7/10/2023)

* Changes to `findAssociatedSeqs`:
    * `groups` argument still exists but is now deprecated and no longer used. Group labels are now automatically determined from the unique values of `group_ids`
    * `sample_ids` argument still exists but is now deprecated and no longer used. Custom sample IDs play no role in `findAssociatedSeqs`; the argument was inherited from a previous function that included the functionality of both `findAssociatedSeqs` and `findAssociatedClones`
* `findPublicClusters` now ignores `plots = TRUE` when `print_plots = FALSE` and `output_dir_unfiltered = NULL`. This prevents unused plots from being generated
* `buildAssociatedClusterNetwork` now uses group ID as the default variable for node colors
* `buildPublicClusterNetwork` and `buildPublicClusterNetworkByRepresentative` now use sample ID as the default variable for node colors
* `buildPublicClusterNetworkByRepresentative` default plot title and subtitle updated for better clarity
* `buildRepSeqNetwork`, `generateNetworkObjects` and `generateNetworkGraphPlots` now use `count_col` as the default variable for node colors if available, followed in priority by cluster ID, then network degree.
* New arguments to `addClusterLabels`:
    * `cluster_id_col` added to permit use with node data where the cluster ID variable has a custom name (e.g., with the output of `buildPublicClusterNetwork`)
    * `greatest_values` added, which can be set to `FALSE` to prioritize the clusters to label based on the least values of the `criterion` variable rather than the greatest values
* A function `exclusiveNodeStats` has been added. This function behaves in the same manner as `chooseNodeStats`, but all arguments are set to `FALSE` by default. Useful when the user only wishes to specify a small number of node-level properties to compute, with all other properties excluded.
* Major revisions to the following vignettes:
    * `NAIR: Network Analysis of Immune Repertoire`
    * `Searching for Public TCR/BCR Clusters`
    * `Searching for Associated TCR/BCR Clusters`
    * `buildRepSeqNetwork`
    * `Network Visualization` (incomplete, in progress)
    * Of particular note, the associated clusters and public clusters vignettes now simulate more reasonable toy data for demonstration purposes.
* `Downstream Analysis` vignette title renamed to `Utility Functions`. A revision to this vignette is planned prior to version 1.0.
* CXX_STD = CXX11 flag removed from src/Makevars and src/Makevars.win
    
    
# Changes in NAIR version 0.0.9036

* `buildRepSeqNetwork` no longer returns an error with `dist_cutoff = 0` (fixed a bug involving the argument checks added in version 0.0.9035).

# 0.0.9035

* Argument checks added to `buildRepSeqNetwork`
* `buildRepSeqNetwork` now automatically attempts to perform the following conversions:
    * coerces the input data to a data frame
    * coerces the sequence column to character
    * coerces the count column to numeric, if provided
* `filterInputData`, which affect top-level functions such as `buildRepSeqNetwork` that call it:
    * automatically drops data rows with `NA` values in the sequence column, with a warning produced
    * added optional `count_col` arg; if provided, the count column will be coerced to numeric and rows with `NA/NaN` values in the count column will be dropped with a warning
* Changes related to choosing node-level network statistics:
    * The `node_stat_settings` function now has a duplicate with the less-confusing name `chooseNodeStats`; the newer name is now used in place of `node_stat_settings` for defaults and in the tutorials
    * The `stats_to_include` argument of `addNodeNetworkStats`, `buildRepSetNetwork`, etc., now also accepts a named logical vector with the same named elements as the list previously required. A list will still work, for backwards compatibility.
    * `chooseNodeStats` / `node_stat_settings` now generate a named logical vector rather than a list.
* The `dist_type` argument is now more flexible in the values it will accept; for example `"lev"` or simply `"l"` is now equivalent to `"levenshtein"`

# Changes in NAIR version 0.0.9034

* Rdocumentation files
    * All examples now use `simulateToyData` to generate data
    
# Changes in NAIR version 0.0.9033

* Rdocumentation files
    * Added documentation for `simulateToyData`
    * All examples now use `simulateToyData` to generate data

# Changes in NAIR version 0.0.9032

* Rdocumentation files
    * Fixed instances where tildes were erroneously used in place of the `\code{}` environment

# Changes in NAIR version 0.0.9031 (2/20/2023)

* Added documentation for the following functions:
    * `addGraphLabels`
    * `addClusterLabels`
* Added the following content to the package vignettes:
    * New vignette `Dual-Chain Network Analysis` for dual-chain network analysis on single-cell data
    * `buildRepSeqNetwork` vignette:
        * `cluster_fun` argument (clustering algorithm)
    * `Network Visualization` vignette:
        * `addClusterLabels` function
    * `Downstream Analysis` vignette:
        * `addClusterLabels` function
    * `Finding Associated Clones` vignette:
        * `addClusterLabels` function used to label the clusters

# Changes in NAIR version 0.0.9030 (1/29/2023)

* Added function `addGraphLabels` for adding text labels to the nodes of a graph plot
* Added function `addClusterLabels` for adding labels to certain clusters in a graph plot

# Changes in NAIR version 0.0.9029 (1/16/2023)

* The algorithm used to identify clusters in `addClusterMembership()` can now be controlled via a new argument `fun`.
* The following functions have a new argument `cluster_fun` that is passed to the `fun` argument of `addClusterMembership()`:
    * `addNodeNetworkStats()`
    * `getClusterStats()`
    * `buildRepSeqNetwork()`
    * `buildAssociatedClusterNetwork()`
    * `findPublicClusters()`
    * `buildPublicClusterNetwork()`
    * `buildPublicClusterNetworkByRepresentative()`

# Changes in NAIR version 0.0.9028 (1/9/2023)

* `buildRepSeqNetwork()` and `generateNetworkObjects()` now return `NULL` with a warning when the constructed network contains no edges.

# Changes in NAIR version 0.0.9027 (1/8/2023)

* `getClusterStats()` now computes sequence-based statistics (e.g., sequence with max count) for dual-chain networks, including a separate set of such statistics for each chain. 
* The name of some variables in the output of `getClusterStats()` have been changed to reflect broader applicability to single-cell data:
    * `max_clone_count` changed to `max_count`
    * `agg_clone_count` changed to `agg_count`

# Changes in NAIR version 0.0.9026 (12/27/2022)

* Added an argument `verbose` to `findAssociatedClones()` that can be optionally set to `TRUE` in order to print additional console output reporting the number of clones in each neighborhood, both by sample and in total.
* Discovered and fixed the following bugs that were present from 0.0.9018 onward:
    * Fixed a bug whereby `findAssociatedSeqs()` was not correctly computing the counts used for Fisher's exact test
    * Fixed a bug in `findPublicClones()` involving identification of the top n clusters by node count in each sample: when more than one cluster possessed the nth highest node count, all of these clusters were included in the top n clusters, resulting in more than n clusters identified by this criterion. This has been reverted to the behavior that existed prior to version 0.0.9018, whereby the first n clusters are selected after sorting data rows by descending node count using the `order` function.

# Changes in NAIR version 0.0.9025 (12/19/2022)

* Fixed a bug in `filterInputData()` that was preventing filtering by minimum sequence length
* Removed `BiocManager` from `Suggests` field of DESCRIPTION, since it is no longer used to access demonstration data when building vignettes.

# Changes in NAIR version 0.0.9024 (12/12/2022)

* Converted all package vignettes to use data created with `simulateToyData()`
* Additional detail added to vignettes for associated clones and public clones workflows


# Changes in NAIR version 0.0.9023 (11/26/2022)

* Added user-level function `simulateToyData` for generating example (toy) data, primarily for use in vignettes, examples and tests.
* Converted README to use `simulateToyData`

# Changes in NAIR version 0.0.9022 (11/25/2022)

* Numerous utility functions that were previously internal have been renamed and exported to be available to the user. These include `generateNetworkObjects()`, `generateNetworkGraphPlots()`, `filterInputData()`, `getNeighborhood()`, `loadDataFromFileList()`, `combineSamples()`, `saveNetwork()`, and `saveNetworkPlots()`.
* Remaining documentation added for all user-level functions.
* The package vignette content has now been split across multiple vignettes, with the package vignette serving as an overview and hub linking to the other vignettes.
* Removed some previously user-facing utility functions from the package's exported namespace that were redundant or unnecessary to expose to the user, including `levAdjacencyMatSparse`, `hamAdjacencyMatSparse`, `generateNetworkFromSeqs`, `getSimilarClones` and `filterClonesBySequenceLength`.



# Changes in NAIR version 0.0.9021 (11/21/2022)

* Minor bug fixes to a few functions in `utils.R` that caused errors or warnings in rare cases


# Changes in NAIR version 0.0.9020 (11/20/2022)

* Internal function `.saveNetwork` changed to user-facing function `saveNetwork`, for use in saving output during downstream analysis

# Changes in NAIR version 0.0.9019 (11/20/2022)

* Bug fixes to associated clones functions

# Changes in NAIR version 0.0.9018 (11/20/2022)

* Changes to `buildRepSeqNetwork()` (many of these changes carry over to other functions):
    * Changed arguments and functionality for saving output. A new `output_type` argument can be used to save the output list to a rds or rda file, rather than the default behavior of saving each item in an individual, uncompressed file. Rather than specifying the filename of each item individually, the `output_name` argument accepts a character string to be used as a common prefix for any files saved. All items are now saved, and the `save_all` argument has been removed.
    * Changed name of argument `other_cols` to `subset_cols` to more accurately reflect its current role (for keeping only certain input columns rather than all)
    * Changed name of argument `drop_chars` to `drop_matches` to better imply that it takes regular expressions and character strings
    * Changed names of arguments `plot_width` and `plot_height` to `pdf_width` and `pdf_height` to more clearly indicate that they affect the dimensions of the saved pdf file, but not those of the plot at the `R` object level (`ggplot`) or as it appears in the R plotting window.
    * Added logical argument `plots` which can be used to prevent plots from being generated.
    * Many of the plotting arguments are now passed to `generateNetworkGraphPlots` via elipses (`...`)
    * The returned value is now always a list and always contains the igraph, adjacency matrix and plots in addition to the node data. The `return_all` argument has been removed.
    * When computing cluster-level properties (`cluster_stats = TRUE`), the corresponding data frame in the output list is now named `cluster_data` (previously was `cluster_stats`)
    * Now invisibly returns its output, such that it can be assigned but will not be printed when the function is called without assigning its output.
    * Now invisibly returns `NULL` with a warning when fewer than 2 sequences exist after filtering; previously it returned an error.
* `buildRepSeqNetwork` now supports a dual-chain approach to analyzing single-cell RepSeq data: two cells (nodes) are considered adjacent if and only if they possess similar receptor sequences in *both* of two chains (e.g., alpha chain and beta chain). This is done by supplying a vector with two column references to `seq_col` instead of a single column reference, where the two columns each contain the receptor sequence from a different chain (e.g., CDR3 sequences from alpha and beta chains) and each row corresponds to a unique cell. This functionality can more generally be used to perform network analysis where similarity is based on any two types of sequences instead of one.
* The functions for the public clusters workflow have been redesigned. Where previously the workflow was handled by a single function, `findPublicClusters`, it is now split across multiple functions in a manner that reduces memory usage and increases the flexibility of the workflow.
    * `findPublicClusters` now performs network analysis on each sample individually to search for public clusters
    * `buildPublicClusterNetwork` combines the public clusters across samples and performs network analysis
    * `buildPublicClusterNetworkByRepresentative` can be used to perform network analysis on the combined public clusters using only a single representative clone from each cluster
    * The optional step for K-means clustering based on numeric encoding of TCR sequences is now performed directly using the `kmeansAtchley` function
* The functions for the associated clusters workflow have been redesigned in a manner similar to those for the public cluster workflow, and now use the same input format as the public cluster functions (one file per sample instead of a single data frame containing all samples):
    * `findAssociatedSeqs` searches across samples for associated clone sequences based on sample membership and Fisher's exact test P-value
    * `findAssociatedClones` searches across samples for clones within a neighborhood of each associated clone sequence
    * `buildAssociatedClusterNetwork` combines the neighborhoods and performs network analysis and clustering
    * Building networks for individual associated clusters/neighborhoods is now done directly using the `buildRepSeqNetwork` function on the desired subset of the output from `buildAssociatedClusterNetwork`
    * K-means clustering on numerically encoded TCR sequences is now done directly using the `kmeansAtchley` function
* A new function `generateNetworkGraphPlots()` has been added, which is capable of generating multiple plots with argument usage similar to that used in `buildRepSeqNetwork` (e.g., multiple color-code variables can be supplied, in which case color scheme and color legend title arguments will meaningfully accept either a scalar or vector valued argument)
* Changes to `plotNetworkGraph()`:
    * Some arguments were renamed and reordered to conform with the arguments of `buildRepSeqNetwork()`
    * Now accepts the argument `show_color_legend = "auto"`, which will show the color legend if `color_nodes_by` is a continuous variable or a discrete variable with at most 20 distinct values.
* `getClusterStats()` can now be used with `seq_col = NULL`, as the sequence variable is only used for a small number of statistics; similar to when `count_col = NULL`, the dependent statistics will be `NA` in the returned data frame, but other cluster properties will still be computed.
* `buildRepSeqNetwork()` and other high-level functions that generate a network from sequences now coerce the list of sequences to a character vector if it is not already in this format (e.g., factors).
* `buildRepSeqNetwork()` and other top-level functions now skip automatic plot generation when more than 1 million nodes are present in the network. This is done to avoid a potential error when calling `ggplot` that occurs when the combined nodes and edges exceed its limitations. After the network is generated and returned, the user can still attempt to manually generate the plot using `plotNetworkGraph()`; in this manner, the potential error will not interfere with completion of building the network.

# Changes in NAIR version 0.0.9017

* `buildDualChainNetwork()` function added


# Changes in NAIR version 0.0.9016

* Package vignette added, which includes an introduction to the package and a tutorial of the `buildRepSeqNetwork()` function
* Package readme file updated
* `findPublicClusters()` now supports `.rds` and `.rda` file types; the `csv_files` argument has been replaced with an argument named `file_type`.
* Identified and fixed an error related to `filterClonesBySequenceLength()` that occurs when the input data only has a single column; this was affecting higher-level functions including `buildRepSeqNetwork()`
* Identified and fixed an error in `getAssociatedClusters()` that occurred when `neighborhood_plots = FALSE` and `return_all = TRUE` (the function tried to include output related to the neighborhood plots when none existed).
* `findAssociatedClones()` now returns an informative error when no sequences pass the filter for minimum sample membership.

# Changes in NAIR version 0.0.9015

* Fixed minor bugs related to changes in 0.0.9014

# Changes in NAIR version 0.0.9014

* The following argument names have been changed in functions in which they appear:
    * `clone_col` to `seq_col`
    * `clones` to `seqs`, (except for `embedClonesByAtchleyFactor()`, for which it was changed to `cdr3_AA`)
    * `edge_dist` to `dist_cutoff`
* The following functions have been renamed:
    * `generateNetworkFromClones` to `generateNetworkFromSeqs`
    * `sparseAdjacencyMatFromClones` to `sparseAdjacencyMatFromSeqs`
    * `adjacencyMatAtchleyFromClones` to `adjacencyMatAtchleyFromSeqs`
    * `embedClonesByAtchleyFactor` to `embedTCRSeqsByAtchleyFactor`
* `getSimilarClones()`: changed default value to of `drop_chars` argument to `NULL`
* `buildRepSeqNetwork()` had its usage revised, primarily regarding the arguments related to the input data:
    * The `nucleo_col`, `amino_col` and `clone_seq_type` arguments have been replaced by a single `seq_col` argument; the function no longer requires both nucleotide and amino acid sequences in the data, and no longer distinguishes between the two
    * The `count_col` is now optional
        * By default, graph plots now use fixed node sizes
    * Other column arguments (`freq_col`, `vgene_col`, `cdr3_length`, etc.) have been removed; these columns were not used for anything specific in the pipeline, so it is not necessary for them to each have their own dedicated argument.
    * By default, all columns of the input data are now carried over to the output. If only some columns are desired, they can be specified using the `other_cols` argument.
    * Input columns are no longer renamed in the output data.
    * The option to aggregate counts/frequencies for identical clones has been removed; this can be done as a data preprocessing step, either manually or using the `aggregateIdenticalClones()` function.
    * An argument `print_plots` has been added to allow the option not to print the plot(s) in `R`. The default is `TRUE`, which corresponds to the previous behavior (all plots are printed).
* `findAssociatedClones()`, `getAssociatedClusters()` and `findPublicClusters()` have had their arguments revised according to the changes to `buildRepSeqNetwork()`.
* R documentation files added for the following functions:
    * `plotNetworkGraph`
    * `sparseAdjacencyMatFromSeqs`
    * `adjacencyMatAtchleyFromSeqs`
    * `embedTCRSeqsByAtchleyFactor`


# Changes in NAIR version 0.0.9013

* R Documentation files added for package functions:
    * `aggregateIdenticalClones`
    * `filterClonesBySequenceLength`
    * `getSimilarClones`
    * `generateNetworkFromClones`
    * `generateNetworkFromAdjacencyMat`
    * `addNodeNetworkStats`
    * `node_stat_settings`
    * `addClusterMembership`
    * `getClusterStats`
* Implemented on-install testing of package functions using package `testthat`


# Changes in NAIR version 0.0.9012

* Package name changed to NAIR (Network Analysis for Immune Repertoire)
* Added `findPublicClusters()`
* `buildClustersAroundSelectedClones()` renamed to `getAssociatedClusters()`
* `getPotentialAssociatedClones()` renamed to `findAssociatedClones()`
* `generateAtchleyCorrHeatmap()` renamed to `kmeansAtchley()`
* `levAdjacencyMatSparse()` and `hamAdjacencyMatSparse()` have a new argument `drop_isolated_nodes` that can be set to `FALSE` to keep isolated nodes. This argument has been added to higher-level functions that dispatch calls to these routines.

# Changes in NAIR version 0.0.9011

* Separate function `embedClonesByAtchleyFactor()` created to perform embedding of TCR CDR3 amino acid sequences in Euclidean 30-space based on Atchley factor representation; this was previously done within the function `adjacencyMatAtchleyFromClones()`, but has now been placed in its own function for more general use
* New function `analyzeDiseaseAssociatedClusters()` created, which is used to perform a combined network analysis on the disease-associated clusters generated by `generateDiseaseAssociatedClusters()`
* New function `generateAtchleyCorrHeatmap()` created
* `graphics`, `reshape2`, `gplots`, `viridisLite` and `RColorBrewer` added as package dependencies via the `Imports` directive of the `DESCRIPTION` file


# Changes in NAIR version 0.0.9010

* `computeMetaForCandidateSeqs()` (helper for `findDiseaseMotifsFromMergedSamples()`) redesigned and renamed to `findDiseaseAssociatedClones()`; this function now takes only the merged sample data as its input data, and filters sequences by a set of criteria (number of samples shared by and minimum seq length) to obtain the list of candidates before conducting the Fisher's exact tests; previously the list of candidates was obtained as input data to the function.
* `findDiseaseMotifsFromMergedSamples()` redesigned to use candidate sequence metadata as input (previously it computed this metadata from the merged sample data) and renamed to `generateDiseaseAssociatedClusters()`
* `generateNetworkWithStats()` now automatically prints the ggraph in R when called (previously the user needed to access the variable `graph_plot` contained in the returned list)
* Added package `dplyr` as a dependency via the `Imports` directive of the `DESCRIPTION` file

# Changes in NAIR version 0.0.9009

* `buildNetwork()` renamed to `generateNetworkWithStats()`
* `adjacencyMatrix()` renamed to `sparseAdjacencyMatFromClones()`
* `genNetworkGraph()` renamed to `generateNetworkFromAdjacencyMat()`
* Added general-purpose helper functions:
    * `aggregateCountsByAminoAcidSeq`
    * `filterDataBySequenceLength`
    * `generateNetworkFromClones`
    * `computeNodeNetworkStats`
    * `addClusterMembership`
    * `computeClusterNetworkStats`
    * `plotNetworkGraph`
    * `adjacencyMatAtchleyFromClones`
* Added function `findDiseaseMotifsFromMergedSamples` and helper functions:
    * `computeMetaForCandidateSeqs`
    * `subsetDataNearTargetMotif`
* Added vignette for `generateNetworkWithStats()`
    
    
# Changes in NAIR version 0.0.9008

* New function `hamDistBounded` for computing bounded Hamming distance in C++
    * Supports strings of unequal length; the longer string is effectively truncated to the length of the shorter string, and the difference in length is added to the Hamming distance with the truncated version.  This is equivalent to extending the shorter string to the length of the longer string by appending placeholder characters that differ from their counterparts in the longer string.
* New function `hamAdjacencyMatSparse` for computing Hamming adjacency matrix in C++
* Changes to function `adjacencyMatrix`: 
    * Now supports `dist_type = "hamming"` to use Hamming distance for determining network adjacency
    * Original indices of input sequences that appear in the adjacency matrix are stored in the row names of the output matrix; the sequences themselves are stored in the column names (accessible via `dimnames()`)
    * The file `col_ids.txt` created by the C++ function that computes the adjacency matrix is now deleted by `adjacencyMatrix()` after it has finished its other tasks.  The information in the file is now stored in the row names of the output matrix and so the file is no longer needed.
    

# Changes in NAIR version 0.0.9007

* `.genNetworkGraph()` internal helper for `buildNetwork()` renamed into a public version `genNetworkGraph()` for use by other package functions and by users; moved to a new file `utils.R` that will be used to house shared helper functions used by multiple package functions.
* Added `inst/python/Atchley_factors.csv`, which stores the Atchley factor amino acid embedding used by `BriseisEncoder.py`


# Changes in NAIR version 0.0.9006

* Python module `tensorflow` added to `installPythonModules()` and Config/reticulate field of the DESCRIPTION file. `tensorflow` is required by the Python module `keras`.


# Changes in NAIR version 0.0.9005

* file `zzz.R`  created with .onUnload() directive to unload package dll via call to `library.dynam.unload()` when the package is unloaded


# Changes in NAIR version 0.0.9004

* Reintroduced reticulate package dependency as well as R .onLoad() directives and R functions related to python integration, dependency management and module installation
* Added Python source code for BriseisEncoder and h5 trained encoder file
* Removed $(SHLIB_OPENMP_CXXFLAGS) from PKG_CXXFLAGS and PKG_LIBS in Makevars as OpenMP is not supported for MacOS (still enabled in Makevars.win). Added comments with instructions for Linux users to enable OpenMP if desired, and added a corresponding note to the main Readme file's installation section for Linux users.


# Changes in NAIR version 0.0.9003

* Removed Python source code (previously used to compute adjacency matrix)
* Removed reticulate package dependency (previously used to manage python dependencies)
* Removed R functions involved in Python management and calling Python scripts
* Changed internal helper R function .createAdjacencyMatrix() to public function adjacencyMatrix() and switched its calls to Python scripts for building the matrix into a call to the corresponding new C++ version added in 0.0.9002 (currently only Levenshtein distance is implemented; Hamming distance will be implemented in a future update)
* Changed buildNetwork() function to use adjacencyMatrix() instead of the previous function .createAdjacencyMatrix()
* Removed levAdjacencyMatDense in favor of always using levAdjacencyMatSparse
* Enabled ARMA_64BIT_WORD to accommodate larger matrices when using the C++ function levAdjacencyMatSparse; specifically, #define ARMA_64BIT_WORD=1 was added to the beginning of each .cpp source file, and -DARMA_64BIT_WORD=1 was added to the PKG_CXXFLAGS definition in Makevars and Makevars.win


# Changes in NAIR version 0.0.9002

* Added internal support for Rcpp
* Added C++ routine for bounded Levenshtein distance (levDistBounded)
* Added C++ routine for dense Levenshtein adjacency matrix (levAdjacencyMatDense)
* Added C++ routine for sparse Levenshtein adjacency matrix (levAdjacencyMatSparse)


# Changes in NAIR version 0.0.9001

* Initial (development) version 
