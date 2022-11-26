
# Inputs:
#       RepSeq data for global cluster network
# Do:
#       embed AA seqs in Euclidean 30-space by Atchley factor
#       perform K-means clustering on embedded values
#       compile data frame where each column corresponds to a sample,
#         each row corresponds to a k-means cluster,
#         each value is the total reads in the sample for all tcr in the cluster
#       Create heatmap for correlation matrix of the above data frame

kmeansAtchley <- function(
  data,
  amino_col = "AminoAcidSeq",
  sample_col = "SampleID",
  group_col,
  k = 100,
  pdf_width = 15,
  pdf_height = 15,
  margin_cluster_heatmap = 25,
  margin_corr_heatmap = 15,
  use_viridis = FALSE, # use viridis color palettes for color-blindness robustness
  output_dir = getwd(),
  file_cluster_heatmap = "atchley_kmeans_cluster_relative_size_profiles_by_sample.pdf",
  file_corr_heatmap = "atchley_kmeans_corr_in_cluster_size_profile_between_samples.pdf",
  return_output = FALSE
) {
  cdr3 <- kmeanClusterID <- NULL
  df <- as.data.frame(data[ , c(amino_col, sample_col, group_col)])
  names(df) <- c("cdr3", "sample_id", "subject_group")

  # Aggregate data by unique rows;
  # add variable UniqueCloneCount to count duplicates of each row
  df <- dplyr::summarize(dplyr::group_by_all(df),
                         UniqueCloneCount = length(cdr3))

  # convert underscores to hyphens in sample ID and subject group values
  # (this allows us to recover the sample ID and group after next step)
  df$sample_id <- gsub("_", "-", df$sample_id)
  df$subject_group <- gsub("_", "-", df$subject_group)

  # Convert data to wide form;
  # first column is cdr3 (one row per unique TCR)
  # one add'l column per sample (append subject group to sample id)
  # cell values in add'l columns are corresponding values of UniqueCloneCount
  df <- reshape2::dcast(df, cdr3 ~ sample_id + subject_group,
                        fun.aggregate = sum, value.var = "UniqueCloneCount")
  cat(paste0(nrow(df), " unique TCR sequences present across all samples.\n"))
  if (k > nrow(df)) { stop("the number 'k' of clusters to build exceeds the number of unique TCR sequences; try a smaller value of 'k'") }

  # Embed amino acid seqs in Euclidean 30-space by Atchley factor representation
  embedded_values <-
    encodeTCRSeqsByAtchleyFactor(df$cdr3, contig_ids = rownames(df))[ , -1]
  rownames(embedded_values) <- df$cdr3

  # Perform K-means clustering on embedded values
  cat(paste0("Performing K-means clustering on the embedded values using k = ", k, " clusters..."))
  set.seed(1)
  clusters_kmeans <- stats::kmeans(embedded_values, k, iter.max = 30)
  df$kmeanClusterID <- clusters_kmeans$cluster # add kmeans cluster ID to data
  cat(" Done.\n")
  num_singleton_clusters <- sum(clusters_kmeans$size == 1)
  if (num_singleton_clusters > 0) {
    warning(paste0("K-means resulted in ", num_singleton_clusters,
                   " clusters containing only a single TCR. If the number of single-TCR clusters is high, try a using a smaller value of 'k', which controls the number of clusters (see help file for details).\n"))
  } else {
    cat(paste0("K-means resulted in ", num_singleton_clusters,
               " clusters containing only a single TCR.\n"))
  }

  # Save cluster ID for each seq
  kmeans_cluster_ids <- df[ , c("cdr3", "kmeanClusterID")]
  names(kmeans_cluster_ids) <- c("AminoAcidSeq", "kmeanClusterID")

  # For each sample, aggregate # of TCR seqs by K-means cluster
  # each row of the data now corresponds to a cluster
  cat("Computing each cluster's share of the TCRs in each sample...")
  df <- dplyr::summarise_all(
    dplyr::group_by(dplyr::select(df, -cdr3), kmeanClusterID), sum)

  # Divide each value by total # of TCRs in the sample
  df <- dplyr::mutate_all(dplyr::select(df, -kmeanClusterID), ~ . / sum(.))
  cat(" Done.\n")

  ## GENERATE HEATMAPS ##
  # Extract subject group of samples from column names of data
  sample_subject_group <-
    sapply(strsplit(colnames(df), "_"), function(x) x[[2]])

  # Remove subject group from column names
  names(df) <- sapply(strsplit(colnames(df), "_"), function(x) x[[1]])

  # Color palette for correlation coefficient (2-color gradient)
  if (use_viridis) {
    colors_corr <- viridisLite::cividis(n = 1000, begin = 1, end = 0)
  } else {
    colors_corr <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "RdBu")))(1000)
  }

  # Color palette for subject group
  num_group_levels <- length(unique(sample_subject_group))
  if (use_viridis) {
    colors_subject_group <-
      viridisLite::viridis(n = num_group_levels)
  } else {
    colors_subject_group <-
      RColorBrewer::brewer.pal(
        n = max(3, num_group_levels), # requires 3 levels minimum
        "Set1")
    if (num_group_levels < 3) {
      colors_subject_group <- colors_subject_group[1:num_group_levels]
    }
  }
  names(colors_subject_group) <- levels(as.factor(data[ , group_col]))

  # Generate heatmap for cluster share of TCRs in each sample
  if (!is.null(file_cluster_heatmap) & !is.null(output_dir)) { # write pdf to file
    cat("Generating a heatmap of each cluster's share of the TCRs in each sample...")
    grDevices::pdf(file.path(output_dir, file_cluster_heatmap),
                   width = pdf_width, height = pdf_height)
    gplots::heatmap.2(as.matrix(df), col = colors_corr, trace = "none",
                      margins = c(margin_cluster_heatmap, margin_cluster_heatmap),
                      lhei = c(1, 3), cexRow = 2, cexCol = 2,
                      key.title = NA, key.xlab = "cluster's share of TCRs in sample",
                      key.par = list(cex = 1.2),
                      ColSideColors = colors_subject_group[sample_subject_group])
    graphics::legend(x = "topright",
                     legend = names(colors_subject_group),
                     col = colors_subject_group,
                     lty = 1, lwd = 10, border = FALSE,
                     bty = "n", y.intersp = 1, cex = 1.75, xpd = TRUE)
    grDevices::dev.off()
    cat(" Done.\n")
    cat(paste0("Plot saved to file:\n  ",
               file.path(output_dir, file_cluster_heatmap), "\n"))
  }

  # Heatmap of correlation between sample profiles of TCR share per cluster
  if (!is.null(file_corr_heatmap) & !is.null(output_dir)) { # write pdf to file
    cat("Generating a heatmap of correlation between samples (based on each sample's profile of TCR shares for all k clusters)...")
    grDevices::pdf(file.path(output_dir, file_corr_heatmap),
                   width = pdf_width, height = pdf_height)
    gplots::heatmap.2(stats::cor(df), col = colors_corr, trace = "none",
                      margins = c(margin_corr_heatmap, margin_corr_heatmap),
                      lhei = c(1, 3), cexRow = 2, cexCol = 2,
                      key.title = NA, key.xlab = "correlation coefficient",
                      key.par = list(cex = 1.2),
                      ColSideColors = colors_subject_group[sample_subject_group])
    graphics::legend(x = "topright",
                     legend = names(colors_subject_group),
                     col = colors_subject_group,
                     lty = 1, lwd = 10, border = FALSE,
                     bty = "n", y.intersp = 1, cex = 1.75, xpd = TRUE)
    grDevices::dev.off()
    cat(" Done.\n")
    cat(paste0("Plot saved to file:\n  ",
               file.path(output_dir, file_corr_heatmap), "\n"))
  }

  if (return_output) {
    return(list("kmeans_cluster_ids" = kmeans_cluster_ids,
                "encoded_values" = embedded_values,
                "cluster_TCR_profiles" = as.matrix(df)))
  }
}