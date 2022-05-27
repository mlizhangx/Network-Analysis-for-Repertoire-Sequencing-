
# Inputs:
#       Combined RepSeq Data From Top Disease-Associated Clusters
# Do:
#       embed AA seqs in Euclidean 30-space by Atchley factor
#       perform K-means clustering on embedded values
#       compile data frame where each column corresponds to a sample,
#         each row corresponds to a k-means cluster,
#         each value is the total reads in the sample for all tcr in the cluster
#       Create heatmap for correlation matrix of the above data frame

generateAtchleyCorrHeatmap <- function(
  data, # combined rep-seq data from disease-associated clusters to be included
  clone_col, # column of data containing clone sequences
  sample_id_col, # column of data containing sample ID
  disease_col, # column of data containing disease status
  k = 100, # number of clusters for K-means clustering
  use_viridis = TRUE, # use viridis color palettes for color-blindness robustness
  output_file = NULL # write pdf to file
) {
  df <- as.data.frame(data[ , c(clone_col, sample_id_col, disease_col)])
  names(df) <- c("cdr3", "sample_id", "disease_status")

  # Aggregate data by unique rows;
  # add variable numReads to count duplicates of each row
  df <- dplyr::summarize(dplyr::group_by_all(df), numReads = length(cdr3))

  # Convert data to wide form;
  # first column is cdr3 (one row per unique TCR CDR3 amino acid seq)
  # one add'l column per sample (append disease status to sample id)
  # cell values in add'l columns are corresponding values of numReads
  df <- reshape2::dcast(df, cdr3 ~ sample_id + disease_status,
                        fun.aggregate = sum, value.var = "numReads")
  cat(paste0(nrow(df), " unique TCR CDR3 amino acid sequences present in data.\n"))

  # Embed amino acid seqs in Euclidean 30-space by Atchley factor representation
  cat("Embedding TCR CDR3 amino acid sequences in 30-dimensional Euclidean space based on Atchley factor representation using trained encoder...\n")
  embedded_values <- embedClonesByAtchleyFactor(df$cdr3,
                                                contig_ids = rownames(df))
  # Perform K-means clustering on embedded values
  cat("Performing K-means clustering on the embedded values...\n")
  set.seed(1)
  clusters_kmeans <- stats::kmeans(embedded_values[ , -1], k, iter.max = 30)
  df$kmean_cluster_id <- clusters_kmeans$cluster # add kmeans cluster ID to data
  cat(paste0("K-means resulted in ", sum(clusters_kmeans$size == 1),
             " clusters containing only a single element. If the number of single-element clusters is high, try a using a smaller value of 'k', which controls the number of clusters.\n"))

  # For each sample, aggregate # of reads across TCR seqs by K-means cluster
  # each row of the data now corresponds to a cluster
  cat("Aggregating the total number of TCR sequence reads by sample and K-means cluster...\n")
  df <- dplyr::summarise_all(
    dplyr::group_by(dplyr::select(df, -cdr3), kmean_cluster_id), sum)

  # Convert aggregate # of reads to share of total reads in the sample (column)
  cat("Transforming aggregate reads into share of sample's total reads...\n")
  df <- dplyr::mutate_all(dplyr::select(df, -kmean_cluster_id), ~ . / sum(.))

  ## GENERATE CORRELATION HEATMAP ##
  cat("Computing sample correlation matrix for vector of cluster-specific aggregate reads across all samples...\n")
  # Extract disease status of samples from column names of data
  sample_disease_status <-
    sapply(strsplit(colnames(df), "_"), function(x) x[[2]])

  # Remove disease status from column names
  names(df) <- sapply(strsplit(colnames(df), "_"), function(x) x[[1]])

  # Color palette for correlation coefficient (2-color gradient)
  if (use_viridis) {
    colors_corr <- viridisLite::cividis(n = 1000, begin = 1, end = 0)
  } else {
    colors_corr <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "RdBu")))(1000)
  }
  # Color palette for disease status
  colors_disease <-
    viridisLite::viridis(n = length(unique(sample_disease_status)))
  names(colors_disease) <- levels(as.factor(data[ , disease_col]))

  if (!is.null(output_file)) { # write pdf to file
    grDevices::pdf(output_file)
    gplots::heatmap.2(stats::cor(df), col = colors_corr, trace = "none",
                      margins = c(15, 15), lhei = c(1, 3),
                      cexRow = 2, cexCol = 2,
                      key.title = NA, key.xlab = "correlation coefficient",
                      key.par = list(cex = 1.2),
                      ColSideColors = colors_disease[sample_disease_status])
    graphics::legend(0.7, 1.125, legend = names(colors_disease),
                     col = colors_disease, lty = 1, lwd = 10, border = FALSE,
                     bty = "n", y.intersp = 0.8, cex = 2, xpd = TRUE)
    grDevices::dev.off()
    cat(paste0("Correlation heatmap saved to file ", output_file, "\n"))
  }

  gplots::heatmap.2(stats::cor(df), col = colors_corr, trace = "none",
                    margins = c(15, 15), lhei = c(1, 3),
                    cexRow = 2, cexCol = 2,
                    key.title = NA, key.xlab = "correlation coefficient",
                    key.par = list(cex = 1.2),
                    ColSideColors = colors_disease[sample_disease_status])
  graphics::legend(0.7, 1.125, legend = names(colors_disease),
                   col = colors_disease, lty = 1, lwd = 10, border = FALSE,
                   bty = "n", y.intersp = 0.8, cex = 2, xpd = TRUE)

  cat("Correlation heatmap printed to plotting window. All tasks complete.\n")
}