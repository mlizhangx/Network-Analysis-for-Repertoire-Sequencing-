library(NAIR)

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


test_that("sparseAdjacencyMatFromSeqs behaves as expected", {

  cloneSeq <- c("A", "AA", "AB", "BB")
  cloneSeqB <- c("A", "AA", "ACCC", "AB", "BB")
  adjMat_k0 <- sparseAdjacencyMatFromSeqs(
    cloneSeq, max_dist = 0, dist_type = "levenshtein"
  )
  adjMat_k1 <- sparseAdjacencyMatFromSeqs(
    cloneSeq, max_dist = 1, dist_type = "levenshtein"
  )
  adjMat_k2 <- sparseAdjacencyMatFromSeqs(
    cloneSeq, max_dist = 2, dist_type = "levenshtein"
  )

  adjMatB_k0 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 0, dist_type = "levenshtein"
  )
  adjMatB_k1 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 1, dist_type = "levenshtein"
  )
  adjMatB_k2 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 2, dist_type = "levenshtein"
  )
  adjMatB_k3 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 3, dist_type = "levenshtein"
  )
  adjMatB_k4 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 4, dist_type = "levenshtein"
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

  adjMat_k0 <- sparseAdjacencyMatFromSeqs(
    cloneSeq, max_dist = 0, dist_type = "hamming"
  )
  adjMat_k1 <- sparseAdjacencyMatFromSeqs(
    cloneSeq, max_dist = 1, dist_type = "hamming"
  )
  adjMat_k2 <- sparseAdjacencyMatFromSeqs(
    cloneSeq, max_dist = 2, dist_type = "hamming"
  )

  adjMatB_k0 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 0, dist_type = "hamming"
  )
  adjMatB_k1 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 1, dist_type = "hamming"
  )
  adjMatB_k2 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 2, dist_type = "hamming"
  )
  adjMatB_k3 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 3, dist_type = "hamming"
  )
  adjMatB_k4 <- sparseAdjacencyMatFromSeqs(
    cloneSeqB, max_dist = 4, dist_type = "hamming"
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
  mat2 <- sparseAdjacencyMatFromSeqs(
    c("fee", "fie", "foe", "fum", "foo")
  )
  expect_s4_class(mat2, "sparseMatrix")
  expect_equal(mat, mat2)
  expect_equal(colnames(mat), colnames(mat2))
  expect_equal(rownames(mat), rownames(mat2))

  mat <- sparseAdjacencyMatFromSeqs(
    c("foo", "foobar", "fubar", "bar")
  )
  expect_s4_class(mat, "sparseMatrix")
  expect_equal(dim(mat), c(0, 0))
  mat <- sparseAdjacencyMatFromSeqs(
    c("foo", "foobar", "fubar", "bar"), max_dist = 2
  )
  expect_s4_class(mat, "sparseMatrix")
  expect_equal(dim(mat), c(0, 0))

  mat <- diag(4)
  mat2 <- sparseAdjacencyMatFromSeqs(
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
  mat2 <- sparseAdjacencyMatFromSeqs(
    c("foo", "foobar", "fubar", "bar"),
    dist_type = "levenshtein",
    max_dist = 2
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
  mat2 <- sparseAdjacencyMatFromSeqs(
    c("foo", "foobar", "fubar", "bar"),
    max_dist = 3
  )
  expect_s4_class(mat2, "sparseMatrix")
  expect_equal(mat, mat2)
  expect_equal(colnames(mat), colnames(mat2))
  expect_equal(rownames(mat), rownames(mat2))
})
