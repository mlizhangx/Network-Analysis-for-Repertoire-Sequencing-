
# Inputs:
#       RepSeq data for global cluster network
# Do:
#       embed AA seqs in Euclidean 30-space by Atchley factor
#       perform K-means clustering on embedded values
#       compile data frame where each column corresponds to a sample,
#         each row corresponds to a k-means cluster,
#         each value is the total reads in the sample for all tcr in the cluster
#       Create heatmap for correlation matrix of the above data frame

generateAtchleyCorrHeatmap <- function(
  data,
  amino_col = "AminoAcidSeq",
  sample_col = "SampleID",
  group_col,
  k = 100,
  plot_width = 15,
  plot_height = 15,
  margin_size = 15,
  use_viridis = FALSE, # use viridis color palettes for color-blindness robustness
  output_dir = getwd(),
  outfile = "atchley_corr_heatmap.pdf", # write pdf to file
  return_output = FALSE
) {

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
  # first column is cdr3 (one row per unique TCR CDR3 amino acid seq)
  # one add'l column per sample (append subject group to sample id)
  # cell values in add'l columns are corresponding values of UniqueCloneCount
  df <- reshape2::dcast(df, cdr3 ~ sample_id + subject_group,
                        fun.aggregate = sum, value.var = "UniqueCloneCount")
  cat(paste0(nrow(df), " unique TCR CDR3 amino acid sequences present in data.\n"))
  if (k > nrow(df)) { stop("the number 'k' of clusters to build exceeds the number of unique TCR sequences; try a smaller value of 'k'") }

  # Embed amino acid seqs in Euclidean 30-space by Atchley factor representation
  embedded_values <- embedClonesByAtchleyFactor(df$cdr3,
                                                contig_ids = rownames(df))[ , -1]
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
                   " clusters containing only a single element. If the number of single-element clusters is high, try a using a smaller value of 'k', which controls the number of clusters (see help file for details).\n"))
  } else {
    cat(paste0("K-means resulted in ", num_singleton_clusters,
               " clusters containing only a single element.\n"))
  }

  # Save cluster ID for each seq
  kmeans_cluster_ids <- df[ , c("cdr3", "kmeanClusterID")]
  names(kmeans_cluster_ids) <- c("AminoAcidSeq", "kmeanClusterID")

  # For each sample, aggregate # of reads across TCR seqs by K-means cluster
  # each row of the data now corresponds to a cluster
  cat("Aggregating the total number of TCR sequence reads for each sample by K-means cluster ID...")
  df <- dplyr::summarise_all(
    dplyr::group_by(dplyr::select(df, -cdr3), kmeanClusterID), sum)
  cat(" Done.\n")

  # Convert aggregate # of reads to share of total reads in the sample (column)
  cat("Transforming aggregate reads into share of sample's total reads...")
  df <- dplyr::mutate_all(dplyr::select(df, -kmeanClusterID), ~ . / sum(.))
  cat(" Done.\nResulting matrix contains one column for each sample and one row for each K-means cluster, with entries encoding the sample's aggregate number of reads for all TCR sequences in the K-means cluster, expressed as a share of the total number of reads for the sample.\n")

  ## GENERATE CORRELATION HEATMAP ##
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
  colors_subject_group <-
    viridisLite::viridis(n = length(unique(sample_subject_group)))
  names(colors_subject_group) <- levels(as.factor(data[ , group_col]))


  heatmap.2(as.matrix(code.df),
            col = hmcol,
            trace = "none",
            margins = c(25,25), lhei = c(1,3),
            cexRow = 1, cexCol = 2,
            key.title = NA, key.xlab = "Cluster frequency",
            key.par = list(cex=1.2),
            ColSideColors = col.color)
  legend(0.7,1.125,
         legend = names(color.key),
         col = color.key,
         lty= 1, lwd = 10,
         border=FALSE, bty="n", y.intersp = 0.8, cex=2, xpd=TRUE)

  heatmap.2(cor_matrix,col = hmcol,
            trace = "none",
            margins = c(15,15), lhei = c(1,3),
            cexRow = 2, cexCol = 2,
            key.title = NA, key.xlab = "correlation coefficient",
            key.par = list(cex=1.2),
            ColSideColors = col.color)
  legend(0.7,1.125,
         legend = names(color.key),
         col = color.key,
         lty= 1, lwd = 10,
         border=FALSE, bty="n", y.intersp = 0.8, cex=2, xpd=TRUE)


  cat("Generating a sample correlation heatmap for the matrix...")
  gplots::heatmap.2(stats::cor(df), col = colors_corr, trace = "none",
                    margins = c(margin_size, margin_size), lhei = c(1, 3),
                    cexRow = 2, cexCol = 2,
                    key.title = NA, key.xlab = "correlation coefficient",
                    key.par = list(cex = 1.2),
                    ColSideColors = colors_subject_group[sample_subject_group])
  graphics::legend(x = 0.75, y = 1.175,
                   legend = names(colors_subject_group),
                   col = colors_subject_group,
                   lty = 1, lwd = 10, border = FALSE,
                   bty = "n", y.intersp = 1, cex = 1.75, xpd = TRUE)
  cat(" Done.\n")
  if (!is.null(outfile) & !is.null(output_dir)) { # write pdf to file
    grDevices::pdf(file.path(output_dir, outfile),
                   width = plot_width, height = plot_height)
    gplots::heatmap.2(stats::cor(df), col = colors_corr, trace = "none",
                      margins = c(margin_size, margin_size), lhei = c(1, 3),
                      cexRow = 2, cexCol = 2,
                      key.title = NA, key.xlab = "correlation coefficient",
                      key.par = list(cex = 1.2),
                      ColSideColors = colors_subject_group[sample_subject_group])
    graphics::legend(x = "topright",
                     legend = names(colors_subject_group),
                     col = colors_subject_group,
                     lty = 1, lwd = 10, border = FALSE,
                     bty = "n", y.intersp = 1, cex = 1.75, xpd = TRUE)
    grDevices::dev.off()
    cat(paste0("Correlation heatmap saved to file:\n  ",
               file.path(output_dir, outfile), "\n"))
  }
  cat("All tasks complete.\n")

  if (return_output) {
    return(list("kmeans_cluster_ids" = kmeans_cluster_ids,
                "embedded_values" = embedded_values))
  }
}