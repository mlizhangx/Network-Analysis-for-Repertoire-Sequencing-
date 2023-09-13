library(NAIR)


# Generate data -----------------------------------------------------------



set.seed(42)
dat <- simulateToyData()
dat[107, c("CloneCount", "CloneFrequency")] <- NA
dat[11, c("CloneCount", "CloneFrequency")] <- NaN
dat[5, c("CloneCount", "CloneFrequency")] <- -Inf
dat[103, c("CloneCount", "CloneFrequency")] <- Inf
dat[191, c("CloneCount", "CloneFrequency")] <- NA
dat[192, c("CloneCount", "CloneFrequency")] <- NaN
dat[193, c("CloneCount", "CloneFrequency")] <- -Inf
dat[198, c("CloneCount", "CloneFrequency")] <- Inf
suppressWarnings(
  net <- buildRepSeqNetwork(dat, "CloneSeq", plots = FALSE)
)
net <- addNodeStats(net)
net <- addClusterStats(net, count_col = "CloneCount",
                       cluster_id_name = "cluster_greedy"
)
net <- addClusterMembership(net = net, cluster_fun = "leiden",
                            cluster_id_name = "cluster_leiden"
)
suppressWarnings(
  net <- addPlots(net, print_plots = FALSE,
                  color_nodes_by = c("cluster_greedy", "cluster_leiden"),
                  color_scheme = "Viridis"
  )
)
suppressWarnings(
  net <- addPlots(net, print_plots = FALSE,
                  color_nodes_by = c("transitivity", "CloneCount"),
                  color_scheme = c("plasma-1", "default"),
                  size_nodes_by = "CloneCount",
                  node_size_limits = c(0.5, 2),
                  color_title = c("Transitivity", "Clone Count"),
                  size_title = "Clone Count"
  )
)
suppressWarnings(
  net <- addPlots(net, print_plots = FALSE,
                  size_nodes_by = 2
  )
)
suppressWarnings(
  net <- addPlots(net, print_plots = FALSE,
                  color_nodes_by = c("eigen_centrality", "authority_score"),
                  size_nodes_by = "eigen_centrality",
                  color_title = c("Eigen Centrality", "auto"),
                  plot_title = NULL,
                  plot_subtitle = NULL
  )
)
suppressWarnings(
  net <- addPlots(net, print_plots = FALSE,
                  color_nodes_by = "coreness",
                  size_nodes_by = "coreness",
                  color_legend = FALSE
  )
)
suppressWarnings(
  net <- addPlots(net, print_plots = FALSE,
                  color_nodes_by = "page_rank",
                  size_nodes_by = "page_rank",
                  color_title = NULL
  )
)
suppressWarnings(
  net2 <- buildRepSeqNetwork(dat, "CloneSeq", print_plots = FALSE,
                             drop_isolated_nodes = FALSE,
                             cluster_stats = TRUE,
                             count_col = "CloneCount"
  )
)
suppressWarnings(
  net3 <- buildRepSeqNetwork(dat, "CloneSeq", print_plots = FALSE,
                             cluster_stats = TRUE
  )
)
suppressWarnings(
  net4 <- buildRepSeqNetwork(dat, "CloneSeq", print_plots = FALSE,
                             node_stats = TRUE
  )
)
suppressWarnings(
  net5 <- buildRepSeqNetwork(dat, "CloneSeq", print_plots = FALSE)
)
dat2 <- simulateToyData(chains = 2)
suppressWarnings(
  sc_net <- buildRepSeqNetwork(
    dat2,
    seq_col = c("AlphaSeq", "BetaSeq"),
    count_col = "UMIs",
    node_stats = TRUE,
    stats_to_include = "all",
    cluster_stats = TRUE,
    color_nodes_by = "SampleID",
    size_nodes_by = "UMIs",
    node_size_limits = c(0.5, 3),
    print_plots = FALSE
  )
)
samples <- simulateToyData(sample_size = 10,
                           output_dir = tempdir()
)
sample_1 <- subset(samples, SampleID == "Sample1")
sample_2 <- subset(samples, SampleID == "Sample2")
utils::write.csv(sample_1, file.path(tempdir(), "Sample1.csv"),
                 row.names = FALSE
)
utils::write.csv(sample_2, file.path(tempdir(), "Sample2.csv"),
                 row.names = FALSE
)
utils::write.csv2(sample_1, file.path(tempdir(), "Sample1b.csv"),
                  row.names = FALSE
)
utils::write.csv2(sample_2, file.path(tempdir(), "Sample2b.csv"),
                  row.names = FALSE
)
utils::write.csv(sample_1, file.path(tempdir(), "Sample1c.csv"),
                 row.names = TRUE
)
utils::write.csv(sample_2, file.path(tempdir(), "Sample2c.csv"),
                 row.names = TRUE
)
utils::write.table(sample_1, file.path(tempdir(), "Sample1.txt"),
                   row.names = TRUE
)
utils::write.table(sample_2, file.path(tempdir(), "Sample2.txt"),
                   row.names = TRUE
)
utils::write.table(sample_1, file.path(tempdir(), "Sample1.tsv"), sep = "\t",
                   row.names = TRUE
)
utils::write.table(sample_2, file.path(tempdir(), "Sample2.tsv"), sep = "\t",
                   row.names = TRUE
)
utils::write.table(sample_1, file.path(tempdir(), "Sample1"),
                   sep = "\t", row.names = TRUE, dec = ",", na = "NA!"
)
utils::write.table(sample_2, file.path(tempdir(), "Sample2"),
                   sep = "\t", row.names = TRUE, dec = ",", na = "NA!"
)
utils::write.table(sample_1, file.path(tempdir(), "Sample1b"),
                   sep = "@", row.names = TRUE, col.names = FALSE
)
utils::write.table(sample_2, file.path(tempdir(), "Sample2b"),
                   sep = "@", row.names = TRUE, col.names = FALSE
)
save(sample_1, file = file.path(tempdir(), "Sample1.rda"))
save(sample_2, file = file.path(tempdir(), "Sample2.rda"))
x <- sample_1
save(x, file = file.path(tempdir(), "Sample1b.rda"))
x <- sample_2
save(x, file = file.path(tempdir(), "Sample2b.rda"))
samples_filtered <- filterInputData(samples,
                                    seq_col = "CloneSeq",
                                    min_seq_length = 13,
                                    drop_matches = "GGG",
                                    subset_cols = c("CloneSeq", "SampleID")
)
samples_filtered$SampleID[samples_filtered$SampleID == "Sample1"] <- "id01"
samples_filtered$SampleID[samples_filtered$SampleID == "Sample2"] <- "id02"

# Distance Functions ------------------------------------------------------



test_that("hamDistBounded works as expected", {

  expect_equal(
    hamDistBounded("foo", "foo", 3),
    0
  )
  expect_equal(
    hamDistBounded("foo", "fee", 3),
    2
  )
  expect_equal(
    hamDistBounded("foo", "fie", 3),
    2
  )
  expect_equal(
    hamDistBounded("foo", "foe", 3),
    1
  )
  expect_equal(
    hamDistBounded("foo", "fum", 3),
    2
  )
  expect_equal(
    hamDistBounded("foo", "bar", 3),
    3
  )
  expect_equal(
    hamDistBounded("foo", "fee", 1),
    -1
  )
  expect_equal(
    hamDistBounded("foo", "fie", 1),
    -1
  )
  expect_equal(
    hamDistBounded("foo", "foe", 1),
    1
  )
  expect_equal(
    hamDistBounded("foo", "fum", 1),
    -1
  )
  expect_equal(
    hamDistBounded("foo", "bar", 1),
    -1
  )
  expect_equal(
    hamDistBounded("foo", "fubar", 10),
    4
  )
  expect_equal(
    hamDistBounded("foo", "foobar", 10),
    3
  )
  expect_equal(
    hamDistBounded("foo", "barfoo", 10),
    6
  )
  expect_equal(
    hamDistBounded("1234567", "1.23457", 7),
    5
  )
  expect_equal(
    hamDistBounded("1234567", "1.23457", 3),
    -1
  )
  expect_equal(
    hamDistBounded("1234567890", "123456789", 10),
    1
  )
  expect_equal(
    hamDistBounded("1234567890", "234567890", 10),
    10
  )
  expect_equal(
    hamDistBounded("foobar", "fubar", 6),
    5
  )
  cloneSeq <- c("A", "AA", "AB", "BB")
  cloneSeqB <- c("A", "AA", "ACCC", "AB", "BB")
  distMat_ham <-
    sapply(cloneSeq,
           function(x) {
             sapply(cloneSeq,
                    function(y) { hamDistBounded(x, y, k = 5) }) })
  distMat_ham_truth <- matrix(1, 4, 4)
  diag(distMat_ham_truth) <- 0
  distMat_ham_truth[c(1, 2), 4] <- 2
  distMat_ham_truth[4, c(1, 2)] <- 2
  expect_equal(distMat_ham, distMat_ham_truth, ignore_attr = TRUE)
  distMatB_ham <-
    sapply(cloneSeqB,
           function(x) {
             sapply(cloneSeqB,
                    function(y) { hamDistBounded(x, y, k = 5) }) })
  distMatB_ham_truth <- matrix(1, 5, 5)
  diag(distMatB_ham_truth) <- 0
  distMatB_ham_truth[c(1, 2), 5] <- 2
  distMatB_ham_truth[5, c(1, 2)] <- 2
  distMatB_ham_truth[3, c(1, 2, 4)] <- 3
  distMatB_ham_truth[c(1, 2, 4), 3] <- 3
  distMatB_ham_truth[5, 3] <- distMatB_ham_truth[3, 5] <- 4
  expect_equal(distMatB_ham, distMatB_ham_truth, ignore_attr = TRUE)


})

test_that("levDistBounded works as expected", {

  expect_equal(
    levDistBounded("foo", "bar", 3),
    3
  )
  expect_equal(
    levDistBounded("1234567", "1.23457", 7),
    2
  )
  expect_equal(
    levDistBounded("1234567", "1.23457", 3),
    2
  )
  expect_equal(
    levDistBounded("1234567890", "123456789", 10),
    1
  )
  expect_equal(
    levDistBounded("1234567890", "234567890", 10),
    1
  )
  expect_equal(
    levDistBounded("foobar", "fubar", 6),
    2
  )
  cloneSeq <- c("A", "AA", "AB", "BB")
  cloneSeqB <- c("A", "AA", "ACCC", "AB", "BB")
  distMat_lev <-
    sapply(cloneSeq,
           function(x) {
             sapply(cloneSeq,
                    function(y) { levDistBounded(x, y, k = 5) }) })
  distMat_lev_truth <- matrix(1, 4, 4)
  diag(distMat_lev_truth) <- 0
  distMat_lev_truth[c(1, 2), 4] <- 2
  distMat_lev_truth[4, c(1, 2)] <- 2
  expect_equal(distMat_lev, distMat_lev_truth, ignore_attr = TRUE)

  # lev distance matrix for cloneSeqB
  distMatB_lev <-
    sapply(cloneSeqB,
           function(x) {
             sapply(cloneSeqB,
                    function(y) { levDistBounded(x, y, k = 5) }) })
  distMatB_lev_truth <- matrix(1, 5, 5)
  diag(distMatB_lev_truth) <- 0
  distMatB_lev_truth[c(1, 2), 5] <- 2
  distMatB_lev_truth[5, c(1, 2)] <- 2
  distMatB_lev_truth[3, c(1, 2, 4)] <- 3
  distMatB_lev_truth[c(1, 2, 4), 3] <- 3
  distMatB_lev_truth[5, 3] <- distMatB_lev_truth[3, 5] <- 4
  expect_equal(distMatB_lev, distMatB_lev_truth, ignore_attr = TRUE)
  expect_equal(levDistBounded("AA", "CACC", 10), 3)
  expect_equal(levDistBounded("CACC", "AA", 10), 3)
  expect_equal(levDistBounded("AA", "ACCC", 10), 3)
  expect_equal(levDistBounded("ACCC", "AA", 10), 3)
  expect_equal(levDistBounded("AA", "CCCA", 10), 3)
  expect_equal(levDistBounded("CCCA", "AA", 10), 3)
  expect_equal(levDistBounded("A", "CCC", 10), 3)
  expect_equal(levDistBounded("CCC", "A", 10), 3)
  expect_equal(levDistBounded("A", "CC", 10), 2)
  expect_equal(levDistBounded("CC", "A", 10), 2)
  expect_equal(10,
               levDistBounded(
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABC*****IJKLMNOPQRST*****Z", 10
               )
  )
  expect_equal(-1,
               levDistBounded(
                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                 "ABC*****IJKLMNOPQRST*****Z", 9
               )
  )
  expect_equal(levDistBounded("", "A", 10), 1)
  expect_equal(levDistBounded("", "", 10), 0)
  expect_equal(levDistBounded("A", "", 10), 1)
  expect_equal(levDistBounded("A", "A", 10), 0)
  expect_equal(levDistBounded("AAAA", "AABA", 10), 1)
  expect_equal(levDistBounded("AAAA", "AABBA", 1), -1)
  expect_equal(levDistBounded("AAAA", "AABBA", 2), 2)
  expect_equal(levDistBounded("AAAA", "", 10), 4)
  expect_equal(levDistBounded("", "AAAA", 10), 4)
})



# Adjacency Matrix --------------------------------------------------------




test_that("generateAdjacencyMatrix behaves as expected", {

  cloneSeq <- c("A", "AA", "AB", "BB")
  cloneSeqB <- c("A", "AA", "ACCC", "AB", "BB")
  expect_warning(
    adjMat_k0 <- generateAdjacencyMatrix(
      cloneSeq, dist_cutoff = 0, dist_type = "levenshtein"
    )
  )
  adjMat_k1 <- generateAdjacencyMatrix(
    cloneSeq, dist_cutoff = 1, dist_type = "levenshtein"
  )
  adjMat_k2 <- generateAdjacencyMatrix(
    cloneSeq, dist_cutoff = 2, dist_type = "levenshtein"
  )

  expect_warning(
    adjMatB_k0 <- generateAdjacencyMatrix(
      cloneSeqB, dist_cutoff = 0, dist_type = "levenshtein"
    )
  )
  adjMatB_k1 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 1, dist_type = "levenshtein"
  )
  adjMatB_k2 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 2, dist_type = "levenshtein"
  )
  adjMatB_k3 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 3, dist_type = "levenshtein"
  )
  adjMatB_k4 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 4, dist_type = "levenshtein"
  )
  adjMat_k1_truth <- adjMat_k2_truth <- adjMatB_k2_truth <-
    matrix(1, nrow = 4, ncol = 4)
  adjMat_k1_truth[c(1, 2), 4] <- 0
  adjMat_k1_truth[4, c(1, 2)] <- 0
  adjMatB_k2_truth <- adjMat_k2_truth
  adjMatB_k1_truth <- adjMat_k1_truth
  adjMatB_k3_truth <- adjMatB_k4_truth <- matrix(1, nrow = 5, ncol = 5)
  adjMatB_k3_truth[3, 5] <- 0
  adjMatB_k3_truth[5, 3] <- 0

  expect_equal(0, length(adjMat_k0))
  expect_equal(
    as.matrix(adjMat_k1), adjMat_k1_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:4, as.numeric(dimnames(adjMat_k1)[[1]])
  )
  expect_equal(
    cloneSeq[1:4], dimnames(adjMat_k1)[[2]]
  )
  expect_equal(
    as.matrix(adjMat_k2), adjMat_k2_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:4, as.numeric(dimnames(adjMat_k2)[[1]])
  )
  expect_equal(
    cloneSeq[1:4], dimnames(adjMat_k2)[[2]]
  )

  expect_equal(0, length(adjMatB_k0))
  expect_equal(
    as.matrix(adjMatB_k1), adjMatB_k1_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    c(1, 2, 4, 5), as.numeric(dimnames(adjMatB_k1)[[1]])
  )
  expect_equal(
    cloneSeqB[c(1, 2, 4, 5)], dimnames(adjMatB_k1)[[2]]
  )
  expect_equal(
    as.matrix(adjMatB_k2), adjMatB_k2_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    c(1, 2, 4, 5), as.numeric(dimnames(adjMatB_k2)[[1]])
  )
  expect_equal(
    cloneSeqB[c(1, 2, 4, 5)], dimnames(adjMatB_k2)[[2]]
  )
  expect_equal(
    as.matrix(adjMatB_k3), adjMatB_k3_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:5, as.numeric(dimnames(adjMatB_k3)[[1]])
  )
  expect_equal(
    cloneSeqB[1:5], dimnames(adjMatB_k3)[[2]]
  )
  expect_equal(
    as.matrix(adjMatB_k4), adjMatB_k4_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:5, as.numeric(dimnames(adjMatB_k4)[[1]])
  )
  expect_equal(
    cloneSeqB[1:5], dimnames(adjMatB_k4)[[2]]
  )

  expect_warning(
    adjMat_k0 <- generateAdjacencyMatrix(
      cloneSeq, dist_cutoff = 0, dist_type = "hamming"
    )
  )
  adjMat_k1 <- generateAdjacencyMatrix(
    cloneSeq, dist_cutoff = 1, dist_type = "hamming"
  )
  adjMat_k2 <- generateAdjacencyMatrix(
    cloneSeq, dist_cutoff = 2, dist_type = "hamming"
  )

  expect_warning(
    adjMatB_k0 <- generateAdjacencyMatrix(
      cloneSeqB, dist_cutoff = 0, dist_type = "hamming"
    )
  )
  adjMatB_k1 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 1, dist_type = "hamming"
  )
  adjMatB_k2 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 2, dist_type = "hamming"
  )
  adjMatB_k3 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 3, dist_type = "hamming"
  )
  adjMatB_k4 <- generateAdjacencyMatrix(
    cloneSeqB, dist_cutoff = 4, dist_type = "hamming"
  )


  expect_equal(0, length(adjMat_k0))
  expect_equal(
    as.matrix(adjMat_k1), adjMat_k1_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:4, as.numeric(dimnames(adjMat_k1)[[1]])
  )
  expect_equal(
    cloneSeq[1:4], dimnames(adjMat_k1)[[2]]
  )
  expect_equal(
    as.matrix(adjMat_k2), adjMat_k2_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:4, as.numeric(dimnames(adjMat_k2)[[1]])
  )
  expect_equal(
    cloneSeq[1:4], dimnames(adjMat_k2)[[2]]
  )

  expect_equal(0, length(adjMatB_k0))
  expect_equal(
    as.matrix(adjMatB_k1), adjMatB_k1_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    c(1, 2, 4, 5), as.numeric(dimnames(adjMatB_k1)[[1]])
  )
  expect_equal(
    cloneSeqB[c(1, 2, 4, 5)], dimnames(adjMatB_k1)[[2]]
  )
  expect_equal(
    as.matrix(adjMatB_k2), adjMatB_k2_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    c(1, 2, 4, 5), as.numeric(dimnames(adjMatB_k2)[[1]])
  )
  expect_equal(
    cloneSeqB[c(1, 2, 4, 5)], dimnames(adjMatB_k2)[[2]]
  )
  expect_equal(
    as.matrix(adjMatB_k3), adjMatB_k3_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:5, as.numeric(dimnames(adjMatB_k3)[[1]])
  )
  expect_equal(
    cloneSeqB[1:5], dimnames(adjMatB_k3)[[2]]
  )
  expect_equal(
    as.matrix(adjMatB_k4), adjMatB_k4_truth,
    ignore_attr = TRUE
  )
  expect_equal(
    1:5, as.numeric(dimnames(adjMatB_k4)[[1]])
  )
  expect_equal(
    cloneSeqB[1:5], dimnames(adjMatB_k4)[[2]]
  )

  mat <- matrix(1, nrow = 4, ncol = 4)
  mat[c(1, 2), 4] <- 0
  mat[4, c(1, 2)] <- 0
  rownames(mat) <- c(1, 2, 3, 5)
  colnames(mat) <- c("fee", "fie", "foe", "foo")
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  mat2 <- generateAdjacencyMatrix(
    c("fee", "fie", "foe", "fum", "foo")
  )
  expect_s4_class(mat2, "sparseMatrix")
  expect_equal(mat, mat2)
  expect_equal(colnames(mat), colnames(mat2))
  expect_equal(rownames(mat), rownames(mat2))

  expect_warning(
    mat <- generateAdjacencyMatrix(
      c("foo", "foobar", "fubar", "bar")
    )
  )
  expect_s4_class(mat, "sparseMatrix")
  expect_equal(dim(mat), c(0, 0))
  expect_warning(
    mat <- generateAdjacencyMatrix(
      c("foo", "foobar", "fubar", "bar"), dist_cutoff = 2
    )
  )
  expect_s4_class(mat, "sparseMatrix")
  expect_equal(dim(mat), c(0, 0))

  mat <- diag(4)
  mat2 <- generateAdjacencyMatrix(
    c("foo", "foobar", "fubar", "bar"),
    drop_isolated_nodes = FALSE
  )
  expect_s4_class(mat2, "sparseMatrix")
  expect_equal(mat, as.matrix(mat2))

  mat <- matrix(1, nrow = 3, ncol = 3)
  mat[1, 3] <- mat[3, 1] <- 0
  rownames(mat) <- c(2, 3, 4)
  colnames(mat) <- c("foobar", "fubar", "bar")
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  mat2 <- generateAdjacencyMatrix(
    c("foo", "foobar", "fubar", "bar"),
    dist_type = "levenshtein",
    dist_cutoff = 2
  )
  expect_s4_class(mat2, "sparseMatrix")
  expect_equal(mat, mat2)
  expect_equal(colnames(mat), colnames(mat2))
  expect_equal(rownames(mat), rownames(mat2))

  mat <- matrix(1, nrow = 3, ncol = 3)
  mat[2, 3] <- mat[3, 2] <- 0
  rownames(mat) <- c(1, 2, 4)
  colnames(mat) <- c("foo", "foobar", "bar")
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  mat2 <- generateAdjacencyMatrix(
    c("foo", "foobar", "fubar", "bar"),
    dist_cutoff = 3
  )
  expect_s4_class(mat2, "sparseMatrix")
  expect_equal(mat, mat2)
  expect_equal(colnames(mat), colnames(mat2))
  expect_equal(rownames(mat), rownames(mat2))

  expect_false(file.exists(file.path(tempdir(), "col_ids.txt")))
})



# Network Building --------------------------------------------------------

test_that("generateNetworkObjects works", {
  expect_true(.isBaseNetworkOutput(generateNetworkObjects(dat, "CloneSeq")))
})

test_that("generateNetworkGraph works", {
  foo <- generateNetworkGraph(net$adjacency_matrix)
  expect_true(.isIgraph(foo))
  expect_true(.doesIgraphMatchData(foo, net$node_data))
})



# Network Analysis --------------------------------------------------------

test_that("addNodeStats works", {
  expect_true(all(names(chooseNodeStats() %in% names(net$node_data))))
  expect_true(.isPosIntegerVector(net$node_data$degree))
  expect_true(.isNumericVector(net$node_data$transitivity))
  expect_true(.hasNAs(net$node_data$transitivity))
  expect_true(.isNumericVector(net$node_data$eigen_centrality))
  expect_true(.isNumericVector(net$node_data$authority_score))
})

test_that("addClusterMembership works", {
  expect_message(
    addClusterMembership(net = net, cluster_id_name = "cluster_greedy"),
    "already contains a variable named"
  )
  expect_true(
    .isPosIntegerVector(net$node_data$cluster_greedy, factor_ok = TRUE)
  )
  expect_true(
    .isPosIntegerVector(net$node_data$cluster_leiden, factor_ok = TRUE)
  )
  expect_equal(
    net$details$clusters_in_network, c("fast_greedy" = 20, "leiden" = 55)
  )
  expect_equal(
    net$details$cluster_id_variable,
    c("fast_greedy" = "cluster_greedy", "leiden" = "cluster_leiden")
  )
})

test_that("addClusterStats works", {
  expect_true(.hasClusterData(net))
  expect_equal(net$details$cluster_data_goes_with, "cluster_greedy")
  expect_equal(net$details$count_col_for_cluster_data, "CloneCount")
  expect_true(all(
    c("cluster_id", "eigen_centrality_eigenvalue") %in% names(net$cluster_data)
  ))
  expect_true(
    .isPosIntegerVector(net$cluster_data$cluster_id, factor_ok = TRUE)
  )
  expect_true(
    .isPosIntegerVector(net$cluster_data$node_count)
  )
  expect_true(
    .isCharVector(net$cluster_data$seq_w_max_count)
  )
  expect_true(
    .isNumericVector(net$cluster_data$diameter_length)
  )
})

test_that("addClusterStats works with NA/Inf values in counts column", {
  expect_equal(net$cluster_data$agg_count[[3]], -Inf)
  expect_equal(net$cluster_data$max_count[[3]], 4422)
  expect_equal(net$cluster_data$agg_count[[6]], Inf)
  expect_equal(net$cluster_data$max_count[[6]], Inf)
  expect_true(.isBaseNetworkOutput(net2))
  expect_true(.hasClusterData(net2))
  expect_equal(net2$cluster_data$agg_count[[3]], -Inf)
  expect_equal(net2$cluster_data$max_count[[3]], 4422)
  expect_equal(net2$cluster_data$agg_count[[6]], Inf)
  expect_equal(net2$cluster_data$max_count[[6]], Inf)
  expect_equal(net2$cluster_data$agg_count[[95]], as.double(NA))
  expect_equal(net2$cluster_data$max_count[[95]], as.double(NA))
  expect_equal(net2$cluster_data$agg_count[[96]], as.double(NA))
  expect_equal(net2$cluster_data$max_count[[96]], as.double(NA))
  expect_equal(net2$cluster_data$agg_count[[97]], -Inf)
  expect_equal(net2$cluster_data$max_count[[97]], -Inf)
  expect_equal(net2$cluster_data$agg_count[[98]], Inf)
  expect_equal(net2$cluster_data$max_count[[98]], Inf)
})

# Loading -----------------------------------------------------------------

test_that("loadDataFromFileList works", {

  expect_error(
    loadDataFromFileList(
      c(tempfile(), tempfile(), tempfile(), tempfile(), tempfile()),
      input_type = "rds"
    ),
    "specifies one or more nonexistent files"
  )

  # RDS files
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, ".rds")),
      input_type = "rds"
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # RDA files, same symbol
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, "b.rda")),
      input_type = "rda",
      data_symbols = "x"
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # RDA files, different symbols
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, ".rda")),
      input_type = "rda",
      data_symbols = c("sample_1", "sample_2")
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # csv files
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, ".csv")),
      input_type = "csv"
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # semicolon-delimited with commas as decimals
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, "b.csv")),
      input_type = "csv2"
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # csv files with row names in first column
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, "c.csv")),
      input_type = "csv",
      read.args = list(row.names = 1)
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # space-separated files
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, ".txt")),
      input_type = "table",
      header = TRUE
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # tab-separated files
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, ".tsv")),
      input_type = "tsv",
      header = TRUE
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # files requiring custom arguments to read.table()
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2)),
      input_type = "table",
      read.args = list(
        header = TRUE,
        sep = "\t",
        dec = ",",
        na.strings = "NA!",
        row.names = 1
      )
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # check that read.args values of header/sep override argument values
  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2)),
      input_type = "table",
      header = FALSE, sep = "",
      read.args = list(
        header = TRUE,
        sep = "\t",
        dec = ",",
        na.strings = "NA!",
        row.names = 1
      )
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  all_samples <-
    loadDataFromFileList(
      file.path(tempdir(), paste0("Sample", 1:2, "b")),
      input_type = "table",
      sep = "@",
      read.args = list(
        row.names = 1,
        col.names = c("rownames",
                      "CloneSeq", "CloneFrequency",
                      "CloneCount", "SampleID"
        )
      )
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)


})

test_that("combineSamples works", {

  # No filters
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".rds")),
      input_type = "rds",
    )
  expect_equal(all_samples, samples, ignore_attr = TRUE)

  # RDS files
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".rds")),
      input_type = "rds",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02"),
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # check subject and group IDs
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".rds")),
      input_type = "rds",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02"),
      subject_ids = c("s01", "s02"),
      group_ids = c("g1", "g2")
    )
  expect_equal(names(all_samples),
               c(names(samples_filtered), "SubjectID", "GroupID"))
  expect_equal(unique(all_samples$SubjectID), c("s01", "s02"))
  expect_equal(unique(all_samples$GroupID), c("g1", "g2"))

  # RDA files, same symbol
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, "b.rda")),
      input_type = "rda",
      data_symbols = "x",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # RDA files, different symbols
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".rda")),
      input_type = "rda",
      data_symbols = c("sample_1", "sample_2"),
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # csv files
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".csv")),
      input_type = "csv",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # semicolon-delimited with commas as decimals
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, "b.csv")),
      input_type = "csv2",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # csv files with row names
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, "c.csv")),
      input_type = "csv",
      read.args = list(row.names = 1),
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # space-separated files
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".txt")),
      input_type = "table",
      header = TRUE,
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # tab-separated files
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".tsv")),
      input_type = "tsv",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2, ".tsv")),
      input_type = "table",
      header = TRUE,
      sep = "\t",
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

  # custom files
  all_samples <-
    combineSamples(
      file.path(tempdir(), paste0("Sample", 1:2)),
      input_type = "table",
      read.args = list(
        header = TRUE, sep = "\t", dec = ",", na.strings = "NA!", row.names = 1
      ),
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGG",
      subset_cols = "CloneSeq",
      sample_ids = c("id01", "id02")
    )
  expect_equal(all_samples, samples_filtered, ignore_attr = TRUE)

})

# Saving ------------------------------------------------------------------

test_that("saveNetwork works", {

  # Individual
  suppressWarnings(saveNetwork(net3, tempdir(), "individual"))
  ndat <- read.csv(file.path(tempdir(), "MyRepSeqNetwork_NodeMetadata.csv"),
                   row.names = 1
  )
  foo <- all.equal(net3$node_data, ndat)
  expect_true(length(foo) == 1)
  expect_match(foo, "is not a factor")
  expect_equal(rownames(net3$node_data), rownames(ndat))
  cdat <- read.csv(file.path(tempdir(), "MyRepSeqNetwork_ClusterMetadata.csv"))
  foo <- all.equal(net3$cluster_data, cdat)
  expect_true(length(foo) == 1)
  expect_match(foo, "is not a factor")
  expect_equal(rownames(net3$cluster_data), rownames(cdat))
  expect_true(
    file.exists(file.path(tempdir(), "MyRepSeqNetwork_Details.txt"))
  )
  expect_true(
    file.exists(file.path(tempdir(), "MyRepSeqNetwork_EdgeList.txt"))
  )
  expect_true(
    file.exists(file.path(tempdir(), "MyRepSeqNetwork_AdjacencyMatrix.mtx"))
  )
  expect_true(
    file.exists(file.path(tempdir(), "MyRepSeqNetwork_Plots.rda"))
  )
  expect_true(
    file.exists(file.path(tempdir(), "MyRepSeqNetwork.pdf"))
  )
  expect_true(
    file.exists(file.path(tempdir(), "MyRepSeqNetwork_GraphLayout.txt"))
  )

  # rds
  suppressWarnings(saveNetwork(net3, tempdir()))
  net3b <- readRDS(file.path(tempdir(), "MyRepSeqNetwork.rds"))
  expect_equal(sapply(net3, typeof), sapply(net3b, typeof))
  expect_true(.isBaseNetworkOutput(net3b))
  expect_true(.hasClusterData(net3b))
  expect_true(.hasPlots(net3b))

  # rda
  suppressWarnings(saveNetwork(net3, tempdir(), "rda"))
  net3_copy <- net3
  load(file.path(tempdir(), "MyRepSeqNetwork.rda"))
  net3b <- net
  net3 <- net3_copy
  expect_equal(sapply(net3, typeof), sapply(net3b, typeof))
  expect_true(.isBaseNetworkOutput(net3b))
  expect_true(.hasClusterData(net3b))
  expect_true(.hasPlots(net3b))
})


test_that("saveNetworkPlots works", {
  suppressWarnings(
    expect_message(
      saveNetworkPlots(net3$plots, outfile = file.path(tempdir(), "plots.pdf"),
                       outfile_layout = file.path(tempdir(), "layout.txt"),
                       verbose = TRUE
      ),
      "saved to file"
    )
  )
  expect_true(file.exists(file.path(tempdir(), "plots.pdf")))
  expect_true(file.exists(file.path(tempdir(), "layout.txt")))
  expect_equal(
    matrix(scan(file.path(tempdir(), "layout.txt"), quiet = TRUE), ncol = 2),
    net3$plots$graph_layout, tolerance = 1e-04
  )

})



# Filtering and Subsetting ------------------------------------------------

test_that("filterInputData works", {
  filtered_data <-
    filterInputData(
      dat,
      seq_col = "CloneSeq",
      min_seq_length = 13,
      drop_matches = "GGGG",
      subset_cols =
        c("CloneSeq", "CloneFrequency", "SampleID")
    )
  expect_equal(nrow(filtered_data), 105)
})

test_that("getNeighborhood works", {
  set.seed(42)
  toy_data <- simulateToyData(sample_size = 500)

  # Get neighborhood around first clone sequence
  nbd <-
    getNeighborhood(
      toy_data,
      seq_col = "CloneSeq",
      target_seq = "GGGGGGGAATTGG"
    )
  expect_true(nrow(nbd) == 67)
  expect_true(ncol(nbd) == 4)
  nbd_t <- table(nbd$CloneSeq)
  expect_true(length(nbd_t) == 22)
  expect_true(max(nbd_t) == 9)
  expect_true(names(nbd_t)[which.max(nbd_t)] == "GGGGGGGAATTGG")
})

test_that("aggregateIdenticalClones works", {
  my_data <- data.frame(
    clone_seq = c("ATCG", rep("ACAC", 2), rep("GGGG", 4)),
    clone_count = rep(1, 7),
    clone_freq = rep(1/7, 7),
    time_point = c("t_0", rep(c("t_0", "t_1"), 3)),
    subject_id = c(rep(1, 5), rep(2, 2))
  )

  data_agg <-
    aggregateIdenticalClones(
      my_data,
      "clone_seq",
      "clone_count",
      "clone_freq",
    )
  expect_equal(dim(data_agg), c(3, 4))

  # group clones by time point
  data_agg_time <-
    aggregateIdenticalClones(
      my_data,
      "clone_seq",
      "clone_count",
      "clone_freq",
      grouping_cols = "time_point"
    )
  expect_equal(dim(data_agg_time), c(5, 5))

  # group clones by subject ID
  data_agg_subject <-
    aggregateIdenticalClones(
      my_data,
      "clone_seq",
      "clone_count",
      "clone_freq",
      grouping_cols = "subject_id"
    )
  expect_equal(dim(data_agg_subject), c(4, 5))

  # group clones by time point and subject ID
  data_agg_time_subject <-
    aggregateIdenticalClones(
      my_data,
      "clone_seq",
      "clone_count",
      "clone_freq",
      grouping_cols =
        c("subject_id", "time_point")
    )
  expect_equal(dim(data_agg_time_subject), c(7, 6))
})



# Visualization -----------------------------------------------------------


test_that("addPlots works", {
  expect_true(.hasPlots(net))
  expect_true(.hasElem(net$plots, "cluster_greedy"))
  expect_true(.hasElem(net$plots, "cluster_leiden"))
  expect_true(.isGgraph(net$plots$cluster_greedy))
  expect_true(.isGgraph(net$plots$cluster_leiden))
})

test_that("extractLayout works", {
  my_layout <- extractLayout(net$plots[[1]])
  expect_equal(my_layout, net$plots$graph_layout, ignore_attr = TRUE)
})







test_that("plots legends behave correctly", {


  expect_equal(names(net$plots$cluster_greedy$guides),
               "colour"
  )
  expect_equal(net$plots$cluster_greedy$guides$colour$title,
               "cluster_greedy"
  )
  expect_equal(net$plots$cluster_greedy$guides$colour$name,
               "legend"
  )

  expect_equal(names(net$plots$cluster_leiden$guides),
               "colour"
  )
  expect_equal(net$plots$cluster_leiden$guides$colour,
               "none"
  )
  expect_match(net$plots$cluster_leiden$labels$subtitle,
               "Nodes colored by cluster_leiden"
  )

  expect_equal(names(net$plots$transitivity$guides),
               c("colour", "size")
  )
  expect_equal(net$plots$transitivity$guides$colour$title,
               "Transitivity"
  )
  expect_equal(net$plots$transitivity$guides$colour$name,
               "colorbar"
  )
  expect_equal(net$plots$transitivity$guides$size$title,
               "Clone Count"
  )
  expect_equal(net$plots$transitivity$guides$size$name,
               "legend"
  )

  expect_equal(names(net$plots$CloneCount$guides),
               c("colour", "size")
  )
  expect_equal(net$plots$CloneCount$guides$colour$title,
               "Clone Count"
  )
  expect_equal(net$plots$CloneCount$guides$colour$name,
               "legend"
  )
  expect_equal(net$plots$CloneCount$guides$size$title,
               "Clone Count"
  )
  expect_equal(net$plots$CloneCount$guides$size$name,
               "legend"
  )

  expect_equal(names(net$plots$eigen_centrality$guides),
               c("colour", "size")
  )
  expect_equal(net$plots$eigen_centrality$guides$colour$title,
               "Eigen Centrality"
  )
  expect_equal(net$plots$eigen_centrality$guides$colour$name,
               "colorbar"
  )
  expect_equal(net$plots$eigen_centrality$guides$size$title,
               "eigen_centrality"
  )
  expect_equal(net$plots$eigen_centrality$guides$size$name,
               "legend"
  )

  expect_equal(names(net$plots$authority_score$guides),
               c("colour", "size")
  )
  expect_equal(net$plots$authority_score$guides$colour$title,
               "authority_score"
  )
  expect_equal(net$plots$authority_score$guides$colour$name,
               "colorbar"
  )
  expect_equal(net$plots$authority_score$guides$size$title,
               "eigen_centrality"
  )
  expect_equal(net$plots$authority_score$guides$size$name,
               "legend"
  )

  expect_equal(names(net$plots$coreness$guides),
               c("colour", "size")
  )
  expect_equal(net$plots$coreness$guides$colour,
               "none"
  )
  expect_equal(net$plots$coreness$guides$size$title,
               "coreness"
  )
  expect_equal(net$plots$coreness$guides$size$name,
               "legend"
  )

  expect_equal(names(net$plots$page_rank$guides),
               c("colour", "size")
  )
  expect_null(net$plots$page_rank$guides$colour$title
  )
  expect_equal(net$plots$page_rank$guides$colour$name,
               "colorbar"
  )
  expect_equal(net$plots$page_rank$guides$size$title,
               "page_rank"
  )
  expect_equal(net$plots$page_rank$guides$size$name,
               "legend"
  )

  expect_equal(names(net2$plots$CloneCount$guides),
               "colour"
  )
  expect_equal(net2$plots$CloneCount$guides$colour$title,
               "CloneCount"
  )
  expect_equal(net2$plots$CloneCount$guides$colour$name,
               "colorbar"
  )

  expect_equal(names(net3$plots), c("cluster_id", "graph_layout"))
  expect_equal(names(net3$plots$cluster_id$guides),
               "colour"
  )
  expect_equal(net3$plots$cluster_id$guides$colour$title,
               "cluster_id"
  )
  expect_equal(net3$plots$cluster_id$guides$colour$name,
               "legend"
  )

  expect_equal(names(net4$plots), c("degree", "graph_layout"))
  expect_equal(names(net4$plots$degree$guides),
               "colour"
  )
  expect_equal(net4$plots$degree$guides$colour$title,
               "degree"
  )
  expect_equal(net4$plots$degree$guides$colour$name,
               "colorbar"
  )

  expect_equal(names(net5$plots), c("uniform_color", "graph_layout"))
  expect_null(names(net5$plots$uniform_color$guides))

  expect_equal(names(sc_net$plots$SampleID$guides),
               c("colour", "size")
  )
  expect_equal(sc_net$plots$SampleID$guides$colour$title,
               "SampleID"
  )
  expect_equal(sc_net$plots$SampleID$guides$colour$name,
               "legend"
  )
  expect_equal(sc_net$plots$SampleID$guides$size$title,
               "UMIs"
  )
  expect_equal(sc_net$plots$SampleID$guides$size$name,
               "legend"
  )

})


test_that("addPlots layout detection works", {
  expect_equal(extractLayout(net$plots$cluster_greedy), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$cluster_leiden), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$cluster_leiden), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$transitivity), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$CloneCount), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$eigen_centrality), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$authority_score), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$coreness), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_equal(extractLayout(net$plots$page_rank), net$plots$graph_layout,
               ignore_attr = TRUE
  )
  expect_match(
    all.equal(extractLayout(net2$plots$CloneCount), net$plots$graph_layout,
              check.attributes = FALSE
    ),
    "400, 244"
  )
  expect_match(
    all.equal(extractLayout(net3$plots[[1]]), net$plots$graph_layout,
              check.attributes = FALSE
    ),
    "1.49"
  )
  expect_match(
    all.equal(extractLayout(net4$plots[[1]]), net$plots$graph_layout,
              check.attributes = FALSE
    ),
    "1.36"
  )
  expect_match(
    all.equal(extractLayout(net5$plots[[1]]), net$plots$graph_layout,
              check.attributes = FALSE
    ),
    "1.49"
  )

  expect_equal(extractLayout(net2$plots[[1]]), net2$plots$graph_layout,
               ignore_attr = TRUE
  )

  expect_equal(extractLayout(net3$plots[[1]]), net3$plots$graph_layout,
               ignore_attr = TRUE
  )

  expect_equal(extractLayout(net4$plots[[1]]), net4$plots$graph_layout,
               ignore_attr = TRUE
  )

  expect_equal(extractLayout(net5$plots[[1]]), net5$plots$graph_layout,
               ignore_attr = TRUE
  )

})


# Dual Chain Analysis -----------------------------------------------------

test_that("buildRepSeqNetwork works with dual-chain", {
  expect_true(.isBaseNetworkOutput(sc_net))
  expect_true(.hasPlots(sc_net))
  expect_true(.hasClusterData(sc_net))
  expect_true(.hasDetails(sc_net))
  expect_equal(sc_net$details$seq_col,
               c(a_col = "AlphaSeq", b_col = "BetaSeq")
  )
  expect_true("mean_A_seq_length" %in% names(sc_net$cluster_data))
})




# Associated Clusters -----------------------------------------------------



set.seed(42)

## Simulate 30 samples from two groups (treatment/control) ##
n_control <- n_treatment <- 15
n_samples <- n_control + n_treatment
sample_size <- 30 # (seqs per sample)
base_seqs <- # first five are associated with treatment
  c("CASSGAYEQYF", "CSVDLGKGNNEQFF", "CASSIEGQLSTDTQYF",
    "CASSEEGQLSTDTQYF", "CASSPEGQLSTDTQYF",
    "RASSLAGNTEAFF", "CASSHRGTDTQYF", "CASDAGVFQPQHF")
# Relative generation probabilities by control/treatment group
pgen_c <- matrix(rep(c(rep(1, 5), rep(30, 3)), times = n_control),
                 nrow = n_control, byrow = TRUE)
pgen_t <- matrix(rep(c(1, 1, rep(1/3, 3), rep(2, 3)), times = n_treatment),
                 nrow = n_treatment, byrow = TRUE)
pgen <- rbind(pgen_c, pgen_t)
simulateToyData(
  samples = n_samples,
  sample_size = sample_size,
  prefix_length = 1,
  prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, n_samples), rep(0, n_samples)),
  affixes = base_seqs,
  affix_probs = pgen,
  num_edits = 0,
  output_dir = tempdir(),
  no_return = TRUE
)
sample_files <-
  file.path(tempdir(),
            paste0("Sample", 1:n_samples, ".rds")
  )
group_labels <- c(rep("reference", n_control),
                  rep("comparison", n_treatment))
associated_seqs <-
  findAssociatedSeqs(
    file_list = sample_files,
    input_type = "rds",
    group_ids = group_labels,
    seq_col = "CloneSeq",
    min_seq_length = NULL,
    drop_matches = NULL,
    min_sample_membership = 0,
    pval_cutoff = 0.1
  )

dir2 <- tempfile()
dir2b <- tempfile()
dir2c <- tempfile()
dir2d <- tempfile()
dir2e <- tempfile()

test_that("findAssociatedSeqs works", {
  expect_true(is.data.frame(associated_seqs))
  expect_true(nrow(associated_seqs) == 4)
  expect_true(ncol(associated_seqs) == 6)

  associated_seqs2 <-
    findAssociatedSeqs(
      file_list = sample_files,
      input_type = "rds",
      group_ids = group_labels,
      subject_ids = rep(1:15, each = 2),
      seq_col = "CloneSeq",
      min_seq_length = NULL,
      drop_matches = NULL,
      min_sample_membership = 0,
      pval_cutoff = 0.1
    )
  expect_true(nrow(associated_seqs2) == 3)
  expect_true(ncol(associated_seqs2) == 9)

})

findAssociatedClones(
  file_list = sample_files,
  input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL,
  drop_matches = NULL,
  output_dir = dir2
)

tmp <- readRDS(list.files(dir2, full.names = TRUE)[[1]])
test_that("findAssociatedClones works", {

  expect_true(length(list.files(dir2)) == 4)
  expect_true(nrow(tmp) == 54)
  expect_true(ncol(tmp) == 6)
})

test_that("buildAssociatedClusterNetwork works", {
  associated_clusters <-
    buildAssociatedClusterNetwork(
      file_list = list.files(dir2,
                             full.names = TRUE
      ),
      seq_col = "CloneSeq",
      size_nodes_by = 1.5,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(associated_clusters))
  expect_true(nrow(associated_clusters$node_data) == 175)
  expect_true(ncol(associated_clusters$node_data) == 16)

})


findAssociatedClones(
  file_list = sample_files,
  input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL,
  drop_matches = NULL,
  output_dir = dir2b,
  output_type = "rda"
)
findAssociatedClones(
  file_list = sample_files,
  input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL,
  drop_matches = NULL,
  output_dir = dir2c,
  output_type = "csv"
)
findAssociatedClones(
  file_list = sample_files,
  input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL,
  drop_matches = NULL,
  output_dir = dir2d,
  output_type = "tsv"
)
findAssociatedClones(
  file_list = sample_files,
  input_type = "rds",
  group_ids = group_labels,
  seq_col = "CloneSeq",
  assoc_seqs = associated_seqs$ReceptorSeq,
  min_seq_length = NULL,
  drop_matches = NULL,
  output_dir = dir2e,
  output_type = "table"
)

test_that("associated clusters functions exchange non-default intermediate data types correctly", {

  associated_clusters <-
    buildAssociatedClusterNetwork(
      file_list = list.files(dir2b,
                             full.names = TRUE
      ),
      input_type = "rda",
      seq_col = "CloneSeq",
      size_nodes_by = 1.5,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(associated_clusters))
  expect_true(nrow(associated_clusters$node_data) == 175)
  expect_true(ncol(associated_clusters$node_data) == 16)


  associated_clusters <-
    buildAssociatedClusterNetwork(
      file_list = list.files(dir2c,
                             full.names = TRUE
      ),
      input_type = "csv",
      seq_col = "CloneSeq",
      size_nodes_by = 1.5,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(associated_clusters))
  expect_true(nrow(associated_clusters$node_data) == 175)
  expect_true(ncol(associated_clusters$node_data) == 16)


  associated_clusters <-
    buildAssociatedClusterNetwork(
      file_list = list.files(dir2d,
                             full.names = TRUE
      ),
      input_type = "tsv",
      seq_col = "CloneSeq",
      size_nodes_by = 1.5,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(associated_clusters))
  expect_true(nrow(associated_clusters$node_data) == 175)
  expect_true(ncol(associated_clusters$node_data) == 16)


  associated_clusters <-
    buildAssociatedClusterNetwork(
      file_list = list.files(dir2e,
                             full.names = TRUE
      ),
      input_type = "table",
      seq_col = "CloneSeq",
      size_nodes_by = 1.5,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(associated_clusters))
  expect_true(nrow(associated_clusters$node_data) == 175)
  expect_true(ncol(associated_clusters$node_data) == 16)

})



# Public Clusters ---------------------------------------------------------



set.seed(42)

## Simulate 30 samples with a mix of public/private sequences ##
n_samples <- 30
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
# Relative generation probabilities
pgen <- cbind(
  stats::toeplitz(0.6^(0:(sample_size - 1))),
  matrix(1, nrow = n_samples, ncol = length(base_seqs) - n_samples)
)
simulateToyData(
  samples = n_samples,
  sample_size = sample_size,
  prefix_length = 1,
  prefix_chars = c("", ""),
  prefix_probs = cbind(rep(1, n_samples), rep(0, n_samples)),
  affixes = base_seqs,
  affix_probs = pgen,
  num_edits = 0,
  output_dir = tempdir(),
  no_return = TRUE
)

sample_files <-
  file.path(tempdir(),
            paste0("Sample", 1:n_samples, ".rds")
  )
dir1 <- tempfile()
findPublicClusters(
  file_list = sample_files,
  input_type = "rds",
  seq_col = "CloneSeq",
  count_col = "CloneCount",
  min_seq_length = NULL,
  drop_matches = NULL,
  top_n_clusters = 3,
  min_node_count = 5,
  min_clone_count = 15000,
  output_dir = dir1
)

tmp <- readRDS(
  list.files(file.path(dir1, "node_meta_data"), full.names = TRUE)[[1]]
)


test_that("buildPublicClusterNetwork works", {
  public_clusters <-
    buildPublicClusterNetwork(
      file_list =
        list.files(
          file.path(dir1, "node_meta_data"),
          full.names = TRUE
        ),
      seq_col = "CloneSeq",
      count_col = "CloneCount",
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(public_clusters))
  expect_true(nrow(public_clusters$node_data) == 517)
  expect_true(ncol(public_clusters$node_data) == 28)
})

test_that("buildPublicClusterNetworkByRepresentative works", {
  public_clusters <-
    buildPublicClusterNetworkByRepresentative(
      file_list =
        list.files(
          file.path(dir1, "cluster_meta_data"),
          full.names = TRUE
        ),
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(public_clusters))
  expect_true(nrow(public_clusters$node_data) == 101)
  expect_true(ncol(public_clusters$node_data) == 31)
})




dir1b <- tempfile()
findPublicClusters(
  file_list = sample_files,
  input_type = "rds",
  seq_col = "CloneSeq",
  count_col = "CloneCount",
  min_seq_length = NULL,
  drop_matches = NULL,
  top_n_clusters = 3,
  min_node_count = 5,
  min_clone_count = 15000,
  output_dir = dir1b,
  output_type = "csv"
)

dir1c <- tempfile()
findPublicClusters(
  file_list = sample_files,
  input_type = "rds",
  seq_col = "CloneSeq",
  count_col = "CloneCount",
  min_seq_length = NULL,
  drop_matches = NULL,
  top_n_clusters = 3,
  min_node_count = 5,
  min_clone_count = 15000,
  output_dir = dir1c,
  output_type = "rda"
)

test_that("public clusters functions exchange non-default intermediate data types correctly", {

  public_clusters <-
    buildPublicClusterNetwork(
      file_list =
        list.files(
          file.path(dir1b, "node_meta_data"),
          full.names = TRUE
        ),
      input_type = "csv",
      seq_col = "CloneSeq",
      count_col = "CloneCount",
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(public_clusters))
  expect_true(nrow(public_clusters$node_data) == 517)
  expect_true(ncol(public_clusters$node_data) == 28)

  public_clusters <-
    buildPublicClusterNetwork(
      file_list =
        list.files(
          file.path(dir1c, "node_meta_data"),
          full.names = TRUE
        ),
      input_type = "rda",
      seq_col = "CloneSeq",
      count_col = "CloneCount",
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(public_clusters))
  expect_true(nrow(public_clusters$node_data) == 517)
  expect_true(ncol(public_clusters$node_data) == 28)


  public_clusters <-
    buildPublicClusterNetworkByRepresentative(
      file_list =
        list.files(
          file.path(dir1c, "cluster_meta_data"),
          full.names = TRUE
        ),
      input_type = "rda",
      size_nodes_by = 1,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(public_clusters))
  expect_true(nrow(public_clusters$node_data) == 101)
  expect_true(ncol(public_clusters$node_data) == 31)

  public_clusters <-
    buildPublicClusterNetworkByRepresentative(
      file_list =
        list.files(
          file.path(dir1b, "cluster_meta_data"),
          full.names = TRUE
        ),
      input_type = "csv",
      size_nodes_by = 1,
      print_plots = FALSE
    )
  expect_true(.isBaseNetworkOutput(public_clusters))
  expect_true(nrow(public_clusters$node_data) == 101)
  expect_true(ncol(public_clusters$node_data) == 31)
})


test_that("findPublicClusters works", {
  expect_true(length(list.files(file.path(dir1, "node_meta_data"))) == 30)
  expect_true(length(list.files(file.path(dir1, "cluster_meta_data"))) == 30)
  # expect_equal(nrow(tmp), 18)  # sometimes reads as 15 (not consistently)
  expect_equal(ncol(tmp), 16)

})

# Clean up temp files -----------------------------------------------------



# clean up temp files
file.remove(
  file.path(
    tempdir(),
    c("MyRepSeqNetwork_NodeMetadata.csv",
      "MyRepSeqNetwork_ClusterMetadata.csv",
      "MyRepSeqNetwork_EdgeList.txt",
      "MyRepSeqNetwork_AdjacencyMatrix.mtx",
      "MyRepSeqNetwork_Details.txt",
      "MyRepSeqNetwork_Plots.rda",
      "MyRepSeqNetwork_GraphLayout.txt",
      "MyRepSeqNetwork.pdf",
      "MyRepSeqNetwork.rds",
      "MyRepSeqNetwork.rda"
    )
  )
)
file.remove(
  file.path(tempdir(),
            c(paste0("Sample", 1:30, ".rds"),
              paste0("Sample", 1:2, ".rda"),
              paste0("Sample", 1:2, "b.rda"),
              paste0("Sample", 1:2, ".csv"),
              paste0("Sample", 1:2, "b.csv"),
              paste0("Sample", 1:2, "c.csv"),
              paste0("Sample", 1:2, ".tsv"),
              paste0("Sample", 1:2, ".txt"),
              paste0("Sample", 1:2),
              paste0("Sample", 1:2, "b")
            )
  )
)
unlink(
  c(dir2, dir2b, dir2c, dir2d, dir2e, dir1, dir1b, dir1c),
  recursive = TRUE
)
