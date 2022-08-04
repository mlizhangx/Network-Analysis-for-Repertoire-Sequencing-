library(NAIR)



# aggregateIdenticalClones ------------------------------------------------

# test_that("aggregateIdenticalClones works with default usage", {
#
#   # Create some data
#   data <- data.frame(
#     clone_seq = c("ATCG", rep("ACAC", 2), rep("GGGG", 4)),
#     clone_count = rep(1, 7),
#     clone_freq = rep(1 / 7, 7),
#     time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
#     subject_id = c(rep(1, 5), rep(2, 2))
#   )
#
#   # Aggregate clones (default usage)
#   data_agg <- aggregateIdenticalClones(
#     data, "clone_seq", "clone_count", "clone_freq")
#
#
#
#   # Aggregate clones by time point
#   data_agg <- aggregateIdenticalClones(
#     data, "clone_seq", "clone_count", "clone_freq",
#     grouping_cols = "time_point")
#
#   # Aggregate clones by subject
#   data_agg <- aggregateIdenticalClones(
#     data, "clone_seq", "clone_count", "clone_freq",
#     grouping_cols = "subject_id")
#
#   # Aggregate clones by subject & time point
#   # (note all clones in each group are already unique)
#   data_agg <- aggregateIdenticalClones(
#     data, "clone_seq", "clone_count", "clone_freq",
#     grouping_cols = c("subject_id", "time_point"))
#
#
#
# })


# # # getSimilarClones
# # Create some data
# data <- data.frame(
#   clone_seq = c("ATCG", "TTCG", "TCG", rep("ATCC", 2), "ATGG", "TGG"),
#   clone_count = rep(1, 7),
#   clone_freq = rep(1 / 7, 7),
#   # group/label variable 1
#   time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
#   # group/label variable 2
#   subject_id = c(rep(1, 5), rep(2, 2))
# )
# target <- "ATCG"
#
# # Default usage (Hamming distance, max distance 1)
# data_target <-
#   getSimilarClones(target, data, "clone_seq")
#
# # Max Hamming distance of 2
# # No sequences have a Hamming distance of 2 from the target sequence,
# # so the output is the same as in the previous example
# data_target2 <-
#   getSimilarClones(target, data, "clone_seq",
#                    max_dist = 2)
#
# # Levenshtein (edit) distance
# data_target_lev <-
#   getSimilarClones(target, data, "clone_seq",
#                    dist_type = "levenshtein")
#
# # Levenshtein (edit) distance, max distance 2
# data_target_lev2 <-
#   getSimilarClones(target, data, "clone_seq",
#                    dist_type = "levenshtein",
#                    max_dist = 2)
#
# # Keeping only the subjects with the target sequence
# data_subjects_lev2 <-
#   getSimilarClones(target, data, "clone_seq",
#                    sample_col = "subject_id",
#                    dist_type = "levenshtein",
#                    max_dist = 2)
#
# # Keeping only the time points with the target sequence
# data_timepoints_lev2 <-
#   getSimilarClones(target, data, "clone_seq",
#                    sample_col = "time_point",
#                    dist_type = "levenshtein",
#                    max_dist = 2)


# # filterClonesBySequenceLength
# # Create some data
# data <- data.frame(
#   clone_seq = c("ATCGATCG", rep("AC", 2), rep("GGGG", 4)),
#   clone_count = rep(1, 7),
#   clone_freq = rep(1 / 7, 7),
#   # group/label variable 1
#   time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
#   # group/label variable 2
#   subject_id = c(rep(1, 5), rep(2, 2))
# )
#
# # Filter using minimum sequence length 3 (default)
# data_filtered <-
#   filterClonesBySequenceLength(data, "clone_seq")
#
# # Filter using minimum sequence length 5
# data_filtered5 <-
#   filterClonesBySequenceLength(data, "clone_seq", min_length = 5)



# generateNetworkFromClones -----------------------------------------------

# # Generate some data
# set.seed(42)
# sample_size <- 200
# group_labels <- rep("Control", times = sample_size)
# assign_group <- sample(c(TRUE, FALSE), size = sample_size, replace = TRUE)
# base_seq_length <- 7
# characters <- c("G", "A", "T", "C")
# char_probs <- c(2/3, 2/9, 1/18, 1/18)
# characters_sample <- sample(characters,
#                             size = base_seq_length * sample_size,
#                             replace = TRUE, prob = char_probs)
# characters_sample <- matrix(characters_sample, ncol = sample_size)
# clone_seqs <- apply(characters_sample, MARGIN = 2,
#                     FUN = function(x) paste0(x, collapse = ""))
# clone_seqs_append <- rep("AATC", times = sample_size)
# append_latent_prob <- runif(sample_size, min = 0, max = 1)
# for (i in 1:sample_size) {
#   case_group <- assign_group[[i]]
#   if (case_group) group_labels[[i]] <- "Case"
#   if ((case_group & append_latent_prob[[i]] > 0.9) |
#       (!case_group & append_latent_prob[[i]] > 0.5)) {
#     clone_seqs_append[[i]] <- "AATCGGGG"
#   } else if ((case_group & append_latent_prob[[i]] > 0.8) |
#              (!case_group & append_latent_prob[[i]] > 0.3)) {
#     clone_seqs_append[[i]] <- "AATCGGT"
#   } else if ((case_group & append_latent_prob[[i]] > 0.7) |
#              (!case_group & append_latent_prob[[i]] > 0.2)) {
#     clone_seqs_append[[i]] <- "AATCGCT"
#   } else if ((case_group & append_latent_prob[[i]] > 0.5) |
#              (!case_group & append_latent_prob[[i]] > 0.1)) {
#     clone_seqs_append[[i]] <- "AATTGCT"
#   } else if ((case_group & append_latent_prob[[i]] > 0.25) |
#              (!case_group & append_latent_prob[[i]] > 0.05)) {
#     clone_seqs_append[[i]] <- "AATTG"
#   }
#   clone_seqs[[i]] <- paste0(clone_seqs[[i]], clone_seqs_append[[i]],
#                             collapse = "")
# }
# counts <- rbinom(sample_size, size = 300, prob = 0.1)
# frequencies <- counts/sum(counts)
# data <- data.frame("clone_seq" = clone_seqs,
#                    "count" = counts,
#                    "frequency" = frequencies,
#                    "group" = group_labels)
#
# net <- generateNetworkFromClones(clones = data$clone_seq,
#                                  dist_type = "levenshtein")
#
#
#
# # Generate some data
# set.seed(42)
# sample_size <- 200
# assign_group <- sample(c(TRUE, FALSE), size = sample_size, replace = TRUE)
# base_seq_length <- 7
# characters <- c("G", "A", "T", "C")
# char_probs <- c(2/3, 2/9, 1/18, 1/18)
# characters_sample <- sample(characters,
#                             size = base_seq_length * sample_size,
#                             replace = TRUE, prob = char_probs)
# characters_sample <- matrix(characters_sample, ncol = sample_size)
# clone_seqs <- apply(characters_sample, MARGIN = 2,
#                     FUN = function(x) paste0(x, collapse = ""))
# clone_seqs_append <- rep("AATC", times = sample_size)
# append_latent_prob <- runif(sample_size, min = 0, max = 1)
# for (i in 1:sample_size) {
#   case_group <- assign_group[[i]]
#   if ((case_group & append_latent_prob[[i]] > 0.9) |
#       (!case_group & append_latent_prob[[i]] > 0.5)) {
#     clone_seqs_append[[i]] <- "AATCGGGG"
#   } else if ((case_group & append_latent_prob[[i]] > 0.8) |
#              (!case_group & append_latent_prob[[i]] > 0.3)) {
#     clone_seqs_append[[i]] <- "AATCGGT"
#   } else if ((case_group & append_latent_prob[[i]] > 0.7) |
#              (!case_group & append_latent_prob[[i]] > 0.2)) {
#     clone_seqs_append[[i]] <- "AATCGCT"
#   } else if ((case_group & append_latent_prob[[i]] > 0.5) |
#              (!case_group & append_latent_prob[[i]] > 0.1)) {
#     clone_seqs_append[[i]] <- "AATTGCT"
#   } else if ((case_group & append_latent_prob[[i]] > 0.25) |
#              (!case_group & append_latent_prob[[i]] > 0.05)) {
#     clone_seqs_append[[i]] <- "AATTG"
#   }
#   clone_seqs[[i]] <- paste0(clone_seqs[[i]], clone_seqs_append[[i]],
#                             collapse = "")
# }
#
# net <- generateNetworkFromClones(clones = clone_seqs,
#                                  dist_type = "levenshtein")
# net <- generateNetworkFromClones(clones = clone_seqs,
#                                  dist_type = "levenshtein",
#                                  return_type = "adjacency_matrix")
# net <- generateNetworkFromClones(clones = clone_seqs,
#                                  dist_type = "euclidean_on_atchley",
#                                  return_type = "adjacency_matrix")
# adjmat <- generateNetworkFromClones(clones = clone_seqs,
#                                     dist_type = "levenshtein",
#                                     return_type = "adjacency_matrix")
# net <- generateNetworkFromAdjacencyMat(adjmat)




# addNodeNetworkStats -----------------------------------------------------

# Generate some data
set.seed(42)
sample_size <- 200
group_labels <- rep("Control", times = sample_size)
assign_group <- sample(c(TRUE, FALSE), size = sample_size, replace = TRUE)
base_seq_length <- 7
characters <- c("G", "A", "T", "C")
char_probs <- c(2/3, 2/9, 1/18, 1/18)
characters_sample <- sample(characters,
                            size = base_seq_length * sample_size,
                            replace = TRUE, prob = char_probs)
characters_sample <- matrix(characters_sample, ncol = sample_size)
clone_seqs <- apply(characters_sample, MARGIN = 2,
                    FUN = function(x) paste0(x, collapse = ""))
clone_seqs_append <- rep("AATC", times = sample_size)
append_latent_prob <- runif(sample_size, min = 0, max = 1)
for (i in 1:sample_size) {
  case_group <- assign_group[[i]]
  if (case_group) group_labels[[i]] <- "Case"
  if ((case_group & append_latent_prob[[i]] > 0.9) |
      (!case_group & append_latent_prob[[i]] > 0.5)) {
    clone_seqs_append[[i]] <- "AATCGGGG"
  } else if ((case_group & append_latent_prob[[i]] > 0.8) |
             (!case_group & append_latent_prob[[i]] > 0.3)) {
    clone_seqs_append[[i]] <- "AATCGGT"
  } else if ((case_group & append_latent_prob[[i]] > 0.7) |
             (!case_group & append_latent_prob[[i]] > 0.2)) {
    clone_seqs_append[[i]] <- "AATCGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.5) |
             (!case_group & append_latent_prob[[i]] > 0.1)) {
    clone_seqs_append[[i]] <- "AATTGCT"
  } else if ((case_group & append_latent_prob[[i]] > 0.25) |
             (!case_group & append_latent_prob[[i]] > 0.05)) {
    clone_seqs_append[[i]] <- "AATTG"
  }
  clone_seqs[[i]] <- paste0(clone_seqs[[i]], clone_seqs_append[[i]],
                            collapse = "")
}
counts <- rbinom(sample_size, size = 300, prob = 0.1)
frequencies <- counts/sum(counts)
data <- data.frame("clone_seq" = clone_seqs,
                   "count" = counts,
                   "frequency" = frequencies,
                   "group" = group_labels)

# Generate network for data
net <- generateNetworkFromClones(data$clone_seq,
                                 drop_isolated_nodes = FALSE)
#
# # Add default network statistics
# data_w_default_stats <- addNodeNetworkStats(data, net)
#
# # Add custom network statistics
# data_w_stats <-
#   addNodeNetworkStats(
#     data, net,
#     stats_to_include =
#       node_stat_settings(
#         cluster_id = TRUE,
#         closeness = TRUE,
#         centrality_by_closeness = TRUE,
#         betweenness = FALSE,
#         centrality_by_betweenness = FALSE,
#         authority_score = FALSE,
#         page_rank = FALSE
#       )
#   )
#
# # Add all network statistics
# data_w_all_stats <-
#   addNodeNetworkStats(
#     data, net,
#     stats_to_include =
#       node_stat_settings(all_stats = TRUE)
#   )
#
#
# # Add cluster ID
# data_w_clusterID <- addClusterMembership(data, net)
#
# # Generate adjacency matrix for data
# adjmat <-
#   generateNetworkFromClones(
#     data$clone_seq,
#     drop_isolated_nodes = FALSE,
#     return_type = "adjacency_matrix"
#   )
#
# # Get cluster stats
# cluster_info <-
#   getClusterStats(data, adjmat,
#                   clone_col = "clone_seq",
#                   count_col = "count")

# Plot network graph
net_plot <- plotNetworkGraph(
  net,
  color_nodes_by = data$group,
  color_scheme = "viridis",
  size_nodes_by = data$count,
  node_size_limits = c(0.5, 5))

# library(AMKAT)
#
# test_that("amkat works with default values", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1)
#   expect_output(str(test1), "List of 15")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_value")
#   expect_match(names(test1)[[12]], "number_of_permutations")
#   expect_match(names(test1)[[13]], "permutation_statistics")
#   expect_match(names(test1)[[14]], "p_value_adjustment")
#   expect_match(names(test1)[[15]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with univariate x", {
#
#   n <- 20; p <- 1; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1)
#   expect_output(str(test1), "List of 14")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "candidate_kernels")
#   expect_match(names(test1)[[9]], "selected_kernels")
#   expect_match(names(test1)[[10]], "test_statistic_value")
#   expect_match(names(test1)[[11]], "number_of_permutations")
#   expect_match(names(test1)[[12]], "permutation_statistics")
#   expect_match(names(test1)[[13]], "p_value_adjustment")
#   expect_match(names(test1)[[14]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 1")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, FALSE)
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with single candidate kernel", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, candidate_kernels = "lin", num_permutations = 1)
#   expect_output(str(test1), "List of 15")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_value")
#   expect_match(names(test1)[[12]], "number_of_permutations")
#   expect_match(names(test1)[[13]], "permutation_statistics")
#   expect_match(names(test1)[[14]], "p_value_adjustment")
#   expect_match(names(test1)[[15]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 1)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with full candidate kernel set", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1,
#                  candidate_kernels = listAmkatKernelFunctions())
#   expect_output(str(test1), "List of 15")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_value")
#   expect_match(names(test1)[[12]], "number_of_permutations")
#   expect_match(names(test1)[[13]], "permutation_statistics")
#   expect_match(names(test1)[[14]], "p_value_adjustment")
#   expect_match(names(test1)[[15]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 5)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
#
# })
# test_that("amkat works without filter", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1, filter_x = FALSE)
#   expect_output(str(test1), "List of 14")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "candidate_kernels")
#   expect_match(names(test1)[[9]], "selected_kernels")
#   expect_match(names(test1)[[10]], "test_statistic_value")
#   expect_match(names(test1)[[11]], "number_of_permutations")
#   expect_match(names(test1)[[12]], "permutation_statistics")
#   expect_match(names(test1)[[13]], "p_value_adjustment")
#   expect_match(names(test1)[[14]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, FALSE)
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
#
# })
# test_that("amkat works with covariates", {
#
#   n <- 20; p <- 2; dim_y <- 3; dim_w <- n - 2
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#   w <- matrix(rnorm(n * dim_w), nrow = n, ncol = dim_w)
#
#   test1 <- amkat(y, x, covariates = w, num_permutations = 1)
#   expect_output(str(test1), "List of 15")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_value")
#   expect_match(names(test1)[[12]], "number_of_permutations")
#   expect_match(names(test1)[[13]], "permutation_statistics")
#   expect_match(names(test1)[[14]], "p_value_adjustment")
#   expect_match(names(test1)[[15]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "int 1")
#   expect_equal(test1$number_of_covariates, 18)
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with floor adjustment to p-value", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1, p_value_adjustment = "floor")
#   expect_output(str(test1), "List of 15")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_value")
#   expect_match(names(test1)[[12]], "number_of_permutations")
#   expect_match(names(test1)[[13]], "permutation_statistics")
#   expect_match(names(test1)[[14]], "p_value_adjustment")
#   expect_match(names(test1)[[15]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment, "Floor value of \\(1 / 1\\) applied")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with no adjustment to p-value", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1, p_value_adjustment = "none")
#   expect_output(str(test1), "List of 15")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_value")
#   expect_match(names(test1)[[12]], "number_of_permutations")
#   expect_match(names(test1)[[13]], "permutation_statistics")
#   expect_match(names(test1)[[14]], "p_value_adjustment")
#   expect_match(names(test1)[[15]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(length(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(test1$selected_kernels), NULL)
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment, "No adjustment")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with multiple test statistics", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1, num_test_statistics = 2)
#   expect_output(str(test1), "List of 17")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "null_residuals")
#   expect_match(names(test1)[[6]], "null_standard_errors")
#   expect_match(names(test1)[[7]], "filter_x")
#   expect_match(names(test1)[[8]], "selected_x_columns")
#   expect_match(names(test1)[[9]], "candidate_kernels")
#   expect_match(names(test1)[[10]], "selected_kernels")
#   expect_match(names(test1)[[11]], "test_statistic_type")
#   expect_match(names(test1)[[12]], "generated_test_statistics")
#   expect_match(names(test1)[[13]], "test_statistic_value")
#   expect_match(names(test1)[[14]], "number_of_permutations")
#   expect_match(names(test1)[[15]], "permutation_statistics")
#   expect_match(names(test1)[[16]], "p_value_adjustment")
#   expect_match(names(test1)[[17]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(is.matrix(test1$null_residuals), TRUE)
#   expect_identical(is.numeric(test1$null_residuals), TRUE)
#   expect_match(typeof(test1$null_residuals), "double")
#   expect_equal(nrow(test1$null_residuals), n)
#   expect_equal(ncol(test1$null_residuals), dim_y)
#   expect_equal(sum(is.na(test1$null_residuals)), 0)
#   expect_equal(sum(is.finite(test1$null_residuals)),
#                length(test1$null_residuals))
#   expect_identical(is.vector(test1$null_standard_errors), TRUE)
#   expect_identical(is.numeric(test1$null_standard_errors), TRUE)
#   expect_match(typeof(test1$null_standard_errors), "double")
#   expect_equal(length(test1$null_standard_errors), dim_y)
#   expect_equal(sum(is.na(test1$null_standard_errors)), 0)
#   expect_equal(sum(is.finite(test1$null_standard_errors)),
#                length(test1$null_standard_errors))
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$selected_x_columns), "double")
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_match(typeof(test1$selected_kernels), "character")
#   expect_equal(nrow(test1$selected_kernels), 2)
#   expect_equal(ncol(test1$selected_kernels), dim_y)
#   expect_identical(.checkCandidateKernels(as.vector(test1$selected_kernels)),
#                    NULL)
#   expect_match(typeof(test1$test_statistic_type), "character")
#   expect_match(typeof(test1$generated_test_statistics), "double")
#   expect_match(typeof(test1$test_statistic_value), "double")
#   expect_equal(sum(is.na(test1$test_statistic_value)), 0)
#   expect_equal(sum(is.finite(test1$test_statistic_value)),
#                length(test1$test_statistic_value))
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(typeof(test1$permutation_statistics), "double")
#   expect_equal(sum(is.na(test1$permutation_statistics)), 0)
#   expect_equal(sum(is.finite(test1$permutation_statistics)),
#                length(test1$permutation_statistics))
#   expect_equal(length(test1$permutation_statistics), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with all individually-optional outputs toggled off", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1, output_test_statistics = FALSE,
#                  output_selected_kernels = FALSE,
#                  output_selected_x_columns = FALSE,
#                  output_null_residuals = FALSE)
#   expect_output(str(test1), "List of 9")
#   expect_match(names(test1)[[1]], "sample_size")
#   expect_match(names(test1)[[2]], "y_dimension")
#   expect_match(names(test1)[[3]], "x_dimension")
#   expect_match(names(test1)[[4]], "number_of_covariates")
#   expect_match(names(test1)[[5]], "filter_x")
#   expect_match(names(test1)[[6]], "candidate_kernels")
#   expect_match(names(test1)[[7]], "number_of_permutations")
#   expect_match(names(test1)[[8]], "p_value_adjustment")
#   expect_match(names(test1)[[9]], "p_value")
#   expect_output(str(test1$sample_size), "int 20")
#   expect_output(str(test1$y_dimension), "int 3")
#   expect_output(str(test1$x_dimension), "int 2")
#   expect_output(str(test1$number_of_covariates), "num 0")
#   expect_identical(test1$filter_x, TRUE)
#   expect_match(typeof(test1$candidate_kernels), "character")
#   expect_equal(length(test1$candidate_kernels), 4)
#   expect_identical(.checkCandidateKernels(test1$candidate_kernels), NULL)
#   expect_equal(length(test1$number_of_permutations), 1)
#   expect_match(test1$p_value_adjustment,
#                "Pseudocount value of \\(1 / 1\\) added")
#   expect_match(typeof(test1$p_value), "double")
#   expect_equal(length(test1$p_value), 1)
#
# })
# test_that("amkat works with output_p_value_only toggled on", {
#
#   n <- 20; p <- 2; dim_y <- 3
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#
#   test1 <- amkat(y, x, num_permutations = 1, output_p_value_only = TRUE)
#   expect_match(typeof(test1), "double")
#   expect_equal(length(test1), 1)
#
# })
# test_that("invalid inputs to amkat are caught and return proper errors", {
#
#   n <- 20; p <- 2; dim_y <- 3;
#   dim_w <- n - 2; dim_w_large <- n - 1
#   n_large <- 25
#   min_sample_size <- 16
#   n_small <- min_sample_size - 1
#   y <- y1 <- y2 <- y3 <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- x1 <- x2 <- x3 <- matrix(rnorm(p * n), nrow = n, ncol = p)
#   y1[1, 1] <- x1[1, 1] <- NA; y2[1, 1] <- x2[1, 1] <- NaN
#   y3[1, 1] <- x3[1, 1] <- Inf; y3[1, 2] <- x3[1, 2] <- -Inf
#   y_small <- matrix(rnorm(dim_y * n_small), nrow = n_small, ncol = dim_y)
#   y_large <- matrix(rnorm(dim_y * n_large), nrow = n_large, ncol = dim_y)
#   w1 <- w2 <- w3 <- matrix(rnorm(n * dim_w), nrow = n, ncol = dim_w)
#   w1[1, 1] <- NA; w2[1, 1] <- NaN; w3[1, 1] <- Inf
#   w_large <- matrix(rnorm(n_large * dim_w), nrow = n_large, ncol = dim_w)
#   w_wide <-  matrix(rnorm(n * dim_w_large), nrow = n, ncol = dim_w_large)
#
#
#   expect_error(amkat(NULL, x), "'y' has zero length")
#   expect_error(amkat(matrix(double()), x), "'y' has zero length")
#   expect_error(amkat(y, NULL), "'x' has zero length")
#   expect_error(amkat(y, matrix(double())), "'x' has zero length")
#   expect_error(amkat(y_small, x)) # nrow(y) < min_sample_size
#   expect_error(amkat(y_large, x),
#                "'y' and 'x' must have the same number of rows")
#   expect_error(amkat(y1, x), "'y' contains NA/NaN values")
#   expect_error(amkat(y2, x), "'y' contains NA/NaN values")
#   expect_error(amkat(y, x1), "'x' contains NA/NaN values")
#   expect_error(amkat(y, x2), "'x' contains NA/NaN values")
#   expect_error(amkat(y3, x), "'y' contains Inf/-Inf values")
#   expect_error(amkat(y, x3), "'x' contains Inf/-Inf values")
#   expect_error(amkat(y, x, covariates = matrix(double())),
#                "'covariates' must either be NULL or have positive length")
#   expect_error(amkat(y, x, covariates = w_large),
#                "'covariates' must have the same row dimension as 'y' and 'x'")
#   expect_error(amkat(y, x, covariates = w1),
#                "'covariates' contains NA/NaN values")
#   expect_error(amkat(y, x, covariates = w2),
#                "'covariates' contains NA/NaN values")
#   expect_error(amkat(y, x, covariates = w3),
#                "'covariates' contains Inf/-Inf values")
#   expect_error(amkat(y, x, covariates = w_wide))
#   expect_error(amkat(y, x, filter_x = logical()),
#                "'filter_x' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(amkat(y, x, filter_x = NULL),
#                "'filter_x' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(amkat(y, x, filter_x = NA),
#                "'filter_x' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(amkat(y, x, filter_x = 42),
#                "'filter_x' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(amkat(y, x, filter_x = c(0, 1)),
#                "'filter_x' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(amkat(y, x, candidate_kernels = character()),
#                "'candidate_kernels' has zero length")
#   expect_error(amkat(y, x, candidate_kernels = NULL),
#                "'candidate_kernels' has zero length")
#   expect_error(amkat(y, x, candidate_kernels = NA),
#                paste0(
#                  "'candidate_kernels' must be a character vector containing ",
#                  "one or more of the following values: \"",
#                  paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))
#   expect_error(amkat(y, x, candidate_kernels = 1),
#                paste0(
#                  "'candidate_kernels' must be a character vector containing ",
#                  "one or more of the following values: \"",
#                  paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))
#   expect_error(amkat(y, x, candidate_kernels = "foo"),
#                paste0(
#                  "'candidate_kernels' must be a character vector containing ",
#                  "one or more of the following values: \"",
#                  paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))
#   expect_error(amkat(y, x, candidate_kernels = c("lin", "g")),
#                paste0(
#                  "'candidate_kernels' must be a character vector containing ",
#                  "one or more of the following values: \"",
#                  paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))
#   expect_error(amkat(y, x, num_permutations = NULL),
#                "'num_permutations' must be a finite, strictly-positive integer")
#   expect_error(amkat(y, x, num_permutations = 0),
#                "'num_permutations' must be a finite, strictly-positive integer")
#   expect_error(amkat(y, x, num_permutations = -1),
#                "'num_permutations' must be a finite, strictly-positive integer")
#   expect_error(amkat(y, x, num_permutations = 1.5),
#                "'num_permutations' must be a finite, strictly-positive integer")
#   expect_error(amkat(y, x, num_permutations = NA),
#                "'num_permutations' must be a finite, strictly-positive integer")
#   expect_error(amkat(y, x, num_permutations = integer()),
#                "'num_permutations' must be a finite, strictly-positive integer")
#   expect_error(amkat(y, x, num_permutations = diag(5)),
#                "'num_permutations' must be a finite, strictly-positive integer")
#
#   expect_error(
#     amkat(y, x, num_test_statistics = NULL),
#     "'num_test_statistics' must be a finite, strictly-positive integer")
#   expect_error(
#     amkat(y, x, num_test_statistics = 0),
#     "'num_test_statistics' must be a finite, strictly-positive integer")
#   expect_error(
#     amkat(y, x, num_test_statistics = -1),
#     "'num_test_statistics' must be a finite, strictly-positive integer")
#   expect_error(
#     amkat(y, x, num_test_statistics = 1.5),
#     "'num_test_statistics' must be a finite, strictly-positive integer")
#   expect_error(
#     amkat(y, x, num_test_statistics = NA),
#     "'num_test_statistics' must be a finite, strictly-positive integer")
#   expect_error(
#     amkat(y, x, num_test_statistics = "foo"),
#     "'num_test_statistics' must be a finite, strictly-positive integer")
#
#   expect_error(
#     amkat(y, x, output_test_statistics = logical()),
#     "'output_test_statistics' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_test_statistics = NULL),
#     "'output_test_statistics' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_test_statistics = NA),
#     "'output_test_statistics' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_test_statistics = 42),
#     "'output_test_statistics' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_test_statistics = c(0, 1)),
#     "'output_test_statistics' must evaluate to \"TRUE\" or \"FALSE\"")
#
#   expect_error(
#     amkat(y, x, output_selected_kernels = logical()),
#     "'output_selected_kernels' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_kernels = NULL),
#     "'output_selected_kernels' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_kernels = NA),
#     "'output_selected_kernels' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_kernels = 42),
#     "'output_selected_kernels' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_kernels = c(0, 1)),
#     "'output_selected_kernels' must evaluate to \"TRUE\" or \"FALSE\"")
#
#   expect_error(
#     amkat(y, x, output_selected_x_columns = logical()),
#     "'output_selected_x_columns' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_x_columns = NULL),
#     "'output_selected_x_columns' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_x_columns = NA),
#     "'output_selected_x_columns' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_x_columns = 42),
#     "'output_selected_x_columns' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_selected_x_columns = c(0, 1)),
#     "'output_selected_x_columns' must evaluate to \"TRUE\" or \"FALSE\"")
#
#   expect_error(
#     amkat(y, x, output_null_residuals = logical()),
#     "'output_null_residuals' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_null_residuals = NULL),
#     "'output_null_residuals' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_null_residuals = NA),
#     "'output_null_residuals' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_null_residuals = 42),
#     "'output_null_residuals' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_null_residuals = c(0, 1)),
#     "'output_null_residuals' must evaluate to \"TRUE\" or \"FALSE\"")
#
#   expect_error(
#     amkat(y, x, output_p_value_only = logical()),
#     "'output_p_value_only' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_p_value_only = NULL),
#     "'output_p_value_only' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_p_value_only = NA),
#     "'output_p_value_only' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_p_value_only = 42),
#     "'output_p_value_only' must evaluate to \"TRUE\" or \"FALSE\"")
#   expect_error(
#     amkat(y, x, output_p_value_only = c(0, 1)),
#     "'output_p_value_only' must evaluate to \"TRUE\" or \"FALSE\"")
#
#   expect_error(amkat(y, x, p_value_adjustment = NULL), paste0(
#     "value of 'p_value_adjustment' must be ",
#     "either \"pseudocount\", \"floor\" or \"none\""))
#   expect_error(amkat(y, x, p_value_adjustment = NA), paste0(
#     "value of 'p_value_adjustment' must be ",
#     "either \"pseudocount\", \"floor\" or \"none\""))
#   expect_error(amkat(y, x, p_value_adjustment = 1), paste0(
#     "value of 'p_value_adjustment' must be ",
#     "either \"pseudocount\", \"floor\" or \"none\""))
#   expect_error(amkat(y, x, p_value_adjustment = 1:5), paste0(
#     "value of 'p_value_adjustment' must be ",
#     "either \"pseudocount\", \"floor\" or \"none\""))
#   expect_error(amkat(y, x, p_value_adjustment = "foo"), paste0(
#     "value of 'p_value_adjustment' must be ",
#     "either \"pseudocount\", \"floor\" or \"none\""))
#   expect_error(amkat(y, x, p_value_adjustment = c("pseudocount", "floor")),
#                paste0("value of 'p_value_adjustment' must be ",
#                       "either \"pseudocount\", \"floor\" or \"none\""))
#
# })