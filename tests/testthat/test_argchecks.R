# library(AMKAT)
#
# # .checkYX ----------------------------------------------------------------------
# test_that(".checkYX properly checks row dim of 'y'", {
#   n <- min_sample_size <- 16; dim_y <- 2; p <- 2
#   n_small <- min_sample_size - 1
#   y <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   y2 <- matrix(rnorm(dim_y * n_small), nrow = n_small, ncol = dim_y)
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#   expect_identical(.checkYX(y, x), NULL)
#   expect_error(.checkYX(y2, x)) # nrow(y) < min_sample_size
# })
#
# test_that(".checkYX properly checks for same row dim of 'x' and 'y'", {
#   n1 <- 20; n2 <- 21; p <- 20; dim_y <- 2
#   y1 <- matrix(rnorm(dim_y * n1), nrow = n1, ncol = dim_y)
#   y2 <- matrix(rnorm(dim_y * n2), nrow = n2, ncol = dim_y)
#   x1 <- matrix(rnorm(p * n1), nrow = n1, ncol = p)
#   x2 <- matrix(rnorm(p * n2), nrow = n2, ncol = p)
#   expect_identical(.checkYX(y1, x1), NULL)
#   expect_identical(.checkYX(y2, x2), NULL)
#   expect_error(.checkYX(y1, x2),
#                "'y' and 'x' must have the same number of rows")
#   expect_error(.checkYX(y2, x1),
#                "'y' and 'x' must have the same number of rows")
# })
#
# test_that(".checkYX properly checks 'x' and 'y' for missing values", {
#   n <- 20; p <- 20; dim_y <- 2
#   y <- y1 <- y2 <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- x1 <- x2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
#   y1[1, 1] <- x1[1, 1] <- NA; y2[1, 1] <- x2[1, 1] <- NaN
#   expect_identical(.checkYX(y, x), NULL)
#   expect_error(.checkYX(y1, x), "'y' contains NA/NaN values")
#   expect_error(.checkYX(y2, x), "'y' contains NA/NaN values")
#   expect_error(.checkYX(y, x1), "'x' contains NA/NaN values")
#   expect_error(.checkYX(y, x2), "'x' contains NA/NaN values")
# })
#
# test_that(".checkYX properly checks 'x' and 'y' for infinite values", {
#   n <- 20; p <- 20; dim_y <- 2
#   y <- y1 <- y2 <- matrix(rnorm(dim_y * n), nrow = n, ncol = dim_y)
#   x <- x1 <- x2 <- matrix(rnorm(p * n), nrow = n, ncol = p)
#   y1[1, 1] <- x1[1, 1] <- Inf; y2[1, 1] <- x2[1, 1] <- -Inf
#   x <- matrix(rnorm(p * n), nrow = n, ncol = p)
#   expect_identical(.checkYX(y, x), NULL)
#   expect_error(.checkYX(y1, x), "'y' contains Inf/-Inf values")
#   expect_error(.checkYX(y2, x), "'y' contains Inf/-Inf values")
#   expect_error(.checkYX(y, x1), "'x' contains Inf/-Inf values")
#   expect_error(.checkYX(y, x2), "'x' contains Inf/-Inf values")
# })
#
# # ------------------------------------------------------------------------------
# test_that("NULL or empty value throws correct error", {
#   expect_identical(.checkNonEmpty("arg", 1), NULL)
#   expect_identical(.checkNonEmpty("arg", ""), NULL)
#   expect_identical(.checkNonEmpty("arg", NA), NULL)
#   expect_error(.checkNonEmpty("arg", NULL), "'arg' has zero length")
#   expect_error(.checkNonEmpty("arg", double()), "'arg' has zero length")
# })
#
# test_that("Non-NULL empty value for 'covariates' throws correct error", {
#   expect_identical(.checkCovariateArgument(1), NULL)
#   expect_identical(.checkCovariateArgument(diag(5)), NULL)
#   expect_identical(.checkCovariateArgument(NULL), NULL)
#   expect_error(.checkCovariateArgument(matrix(double())),
#                "'covariates' must either be NULL or have positive length")
# })
#
# test_that("argument not evaluating to TRUE or FALSE throws correct error", {
#   expectNotBoolean <- function(x) {
#     eval(bquote(
#       expect_error(.checkTrueOrFalse("argument", .(x)),
#                    "'argument' must evaluate to \"TRUE\" or \"FALSE\"")))
#   }
#   expect_identical(.checkTrueOrFalse("arg", TRUE), NULL)
#   expect_identical(.checkTrueOrFalse("arg", FALSE), NULL)
#   expectNotBoolean(logical())
#   expectNotBoolean(NULL)
#   expectNotBoolean(NA)
#   expectNotBoolean(c(NA, TRUE))
#   expectNotBoolean(c(TRUE, FALSE))
#   expectNotBoolean(c(0, 1))
#   expectNotBoolean(42)
# })
#
# test_that("invalid input for 'candidate_kernels' throws correct error", {
#   expect_identical(.checkCandidateKernels(listAmkatKernelFunctions()), NULL)
#   expect_identical(.checkCandidateKernels("lin"), NULL)
#   expect_identical(.checkCandidateKernels("quad"), NULL)
#   expect_identical(.checkCandidateKernels("gau"), NULL)
#   expect_identical(.checkCandidateKernels("exp"), NULL)
#   expect_identical(.checkCandidateKernels("IBS"), NULL)
#   expect_identical(.checkCandidateKernels(c("lin", "IBS")), NULL)
#   expect_identical(
#     .checkCandidateKernels(c("lin", "quad", "gau", "exp", "IBS")), NULL)
#   expect_error(.checkCandidateKernels(character()),
#                "'candidate_kernels' has zero length")
#   expect_error(.checkCandidateKernels(NULL),
#                "'candidate_kernels' has zero length")
#   expectNotKernels <- function(x) {
#     eval(bquote(expect_error(.checkCandidateKernels(.(x)), paste0(
#       "'candidate_kernels' must be a character vector containing ",
#       "one or more of the following values: \"",
#       paste(listAmkatKernelFunctions(), collapse = "\", \""), "\""))))
#   }
#   expectNotKernels(NA)
#   expectNotKernels(1)
#   expectNotKernels(1:5)
#   expectNotKernels("foo")
#   expectNotKernels(c("lin", "foo"))
#   expectNotKernels(c("lin", NA))
# })
#
# test_that("value other than positive integer throws correct error", {
#   expect_identical(.checkPositiveInteger('arg', 1), NULL)
#   expect_identical(.checkPositiveInteger('arg', 2), NULL)
#   expect_identical(.checkPositiveInteger('arg', 42), NULL)
#   expect_identical(.checkPositiveInteger('arg', 9001), NULL)
#   expectNotPositiveInteger <- function(x) {
#     eval(bquote(expect_error(.checkPositiveInteger('arg', .(x)), paste0(
#       "'arg' must be a finite, strictly-positive integer"))))
#   }
#   expectNotPositiveInteger(0)
#   expectNotPositiveInteger(-1)
#   expectNotPositiveInteger(1.5)
#   expectNotPositiveInteger("foo")
#   expectNotPositiveInteger(c(1,2))
#   expectNotPositiveInteger(Inf)
#   expectNotPositiveInteger(NA)
#   expectNotPositiveInteger(NULL)
#   expectNotPositiveInteger(matrix(1, 2, 2))
#   expectNotPositiveInteger(integer())
# })
#
# test_that("invalid input for 'p_value_adjustment' throws correct error", {
#   expect_identical(.checkPValueAdjustment('pseudocount'), NULL)
#   expect_identical(.checkPValueAdjustment('floor'), NULL)
#   expect_identical(.checkPValueAdjustment('none'), NULL)
#   expectNotPValueAdjustment <- function(x) {
#     eval(bquote(expect_error(.checkPValueAdjustment(.(x)), paste0(
#       "value of 'p_value_adjustment' must be ",
#       "either \"pseudocount\", \"floor\" or \"none\""))))
#   }
#   expectNotPValueAdjustment(NULL)
#   expectNotPValueAdjustment(NA)
#   expectNotPValueAdjustment(1)
#   expectNotPValueAdjustment(1:5)
#   expectNotPValueAdjustment(character())
#   expectNotPValueAdjustment("foo")
#   expectNotPValueAdjustment(c("pseudocount", "floor"))
#   expectNotPValueAdjustment(c("pseudocount", NA))
# })
#
# test_that("'covariates' with NA/Inf values or bad dim throws correct error", {
#   w <- matrix(rnorm(2 * 20), nrow = 20, ncol = 2)
#   w0 <- w1 <- w2 <- w3 <- matrix(rnorm(15 * 16), nrow = 16, ncol = 15)
#   w1[1, 1] <- NA; w2[1, 1] <- NaN; w3[1, 1] <- Inf
#   expect_identical(.checkCovariateContent(w, 20), NULL)
#   expect_error(.checkCovariateContent(w0, 20),
#                "'covariates' must have the same row dimension as 'y' and 'x'")
#   expect_error(.checkCovariateContent(w1, 16),
#                "'covariates' contains NA/NaN values")
#   expect_error(.checkCovariateContent(w2, 16),
#                "'covariates' contains NA/NaN values")
#   expect_error(.checkCovariateContent(w3, 16),
#                "'covariates' contains Inf/-Inf values")
#   expect_error(.checkCovariateContent(w0, 16)) # ncol(w0) > n - 2
# })
