library(NAIR)


# Generate data -----------------------------------------------------------



set.seed(42)
dat <- simulateToyData()
suppressWarnings(
  net <- buildRepSeqNetwork(dat, "CloneSeq", cluster_stats = TRUE)
)

suppressWarnings(
  net0 <- buildRepSeqNetwork(dat, "CloneSeq", drop_isolated_nodes = FALSE)
)
net1 <- net2 <- net0
net1$plots$graph_layout <- NULL
net2$plots$graph_layout <- net$plots$graph_layout


# Generic Checks ----------------------------------------------------------


test_that(".orNull works correctly", {
  expect_error(
    .orNull(.noNAs, NA, "argument"),
    "must not contain NA or NaN values"
  )
  expect_error(
    .orNull(.MUST.hasLength1, 1:3, "argument"),
    "must have length 1"
  )
  expect_error(
    .orNull(.MUST.isNumeric, "foo", "argument"),
    paste("must be of type", dQuote("numeric"))
  )
  expect_error(
    .orNull(.MUST.isNumeric, "foo", "argument"),
    paste("must be of type", dQuote("numeric"))
  )
  expect_true(.orNull(.noNAs, NULL))
  expect_true(.orNull(.MUST.hasLength1, NULL))
  expect_true(.orNull(.MUST.isNumeric, NULL))
})


# Universal Properties ----------------------------------------------------

test_that(".MUST.hasLength works correctly", {
  expect_error(
    .MUST.hasLength(NULL, 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(logical(0), 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(numeric(0), 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(character(0), 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(diag(2), 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(1:3, 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(3, 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength("foo", 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(NA, 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(NaN, 2),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength(Inf, 2),
    "must have length 2"
  )
  expect_null(.MUST.hasLength(c("fee", "fie"), 2))
  expect_null(.MUST.hasLength(c(23, NaN), 2))
  expect_null(.MUST.hasLength(c(TRUE, FALSE), 2))
  expect_null(.MUST.hasLength(numeric(2), 2))
  expect_error(
    .MUST.hasLength(NULL, 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(logical(0), 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(numeric(0), 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(character(0), 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(c("fee", "fie"), 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(diag(2), 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(1:3, 1),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength(1:3, c(1, 2)),
    "must have one of the following lengths: 1, 2"
  )
  expect_null(.MUST.hasLength(3, 1))
  expect_null(.MUST.hasLength(3, c(1, 2)))
  expect_null(.MUST.hasLength(3:4, c(1, 2)))
  expect_null(.MUST.hasLength("foo", 1))
  expect_null(.MUST.hasLength(NA, 1))
  expect_null(.MUST.hasLength(NaN, 1))
  expect_null(.MUST.hasLength(Inf, 1))
  expect_error(
    .MUST.hasLength(1:4, 3),
    "must have length 3"
  )
  expect_null(.MUST.hasLength(1:3, 3))
})

test_that(".MUST.hasPosLength works correctly", {
  expect_error(
    .MUST.hasPosLength(NULL, "argument"),
    "must have positive length"
  )
  expect_error(
    .MUST.hasPosLength(logical(0), "argument"),
    "must have positive length"
  )
  expect_error(
    .MUST.hasPosLength(numeric(0), "argument"),
    "must have positive length"
  )
  expect_error(
    .MUST.hasPosLength(character(0), "argument"),
    "must have positive length"
  )
  expect_null(.MUST.hasPosLength(3, "argument"))
  expect_null(.MUST.hasPosLength("foo", "argument"))
  expect_null(.MUST.hasPosLength(NA, "argument"))
  expect_null(.MUST.hasPosLength(NaN, "argument"))
  expect_null(.MUST.hasPosLength(Inf, "argument"))
  expect_null(.MUST.hasPosLength(TRUE, "argument"))
  expect_null(.MUST.hasPosLength(c("fee", "fie"), "argument"))
  expect_null(.MUST.hasPosLength(c(23, NaN), "argument"))
  expect_null(.MUST.hasPosLength(c(TRUE, FALSE), "argument"))
  expect_null(.MUST.hasPosLength(numeric(2), "argument"))
  expect_null(.MUST.hasPosLength(1:4, "argument"))
  expect_null(.MUST.hasPosLength(3, "argument"))
  expect_null(.MUST.hasPosLength(diag(2), "argument"))
})

test_that(".MUST.hasLength1 works correctly", {
  expect_error(
    .MUST.hasLength1(NULL, "argument"),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength1(logical(0), "argument"),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength1(numeric(0), "argument"),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength1(character(0), "argument"),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength1(c("fee", "fie"), "argument"),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength1(diag(2), "argument"),
    "must have length 1"
  )
  expect_error(
    .MUST.hasLength1(1:3, "argument"),
    "must have length 1"
  )
  expect_null(.MUST.hasLength1(3, "argument"))
  expect_null(.MUST.hasLength1("foo", "argument"))
  expect_null(.MUST.hasLength1(NA, "argument"))
  expect_null(.MUST.hasLength1(NaN, "argument"))
  expect_null(.MUST.hasLength1(Inf, "argument"))
  expect_null(.MUST.hasLength1(TRUE, "argument"))
})

test_that(".MUST.hasLength2 works correctly", {
  expect_error(
    .MUST.hasLength2(NULL, "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(logical(0), "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(numeric(0), "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(character(0), "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(diag(2), "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(1:3, "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(3, "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2("foo", "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(NA, "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(NaN, "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(Inf, "argument"),
    "must have length 2"
  )
  expect_error(
    .MUST.hasLength2(TRUE, "argument"),
    "must have length 2"
  )
  expect_null(.MUST.hasLength2(c("fee", "fie"), "argument"))
  expect_null(.MUST.hasLength2(c(23, NaN), "argument"))
  expect_null(.MUST.hasLength2(c(TRUE, FALSE), "argument"))
  expect_null(.MUST.hasLength2(numeric(2), "argument"))
})

test_that(".MUST.hasElement works correctly", {
  foo <- c("first" = 1, "second" = 2, "third" = 3)
  expect_error(
    .MUST.hasElement(foo, "fourth"),
    "must contain an element named"
  )
  expect_error(
    .MUST.hasElement(NULL, "first"),
    "must contain an element named"
  )
  expect_error(
    .MUST.hasElement(1:3, "first"),
    "must contain an element named"
  )
  expect_error(
    .MUST.hasElement(c("first", "second", "third"), "first"),
    "must contain an element named"
  )
  expect_null(.MUST.hasElement(foo, "first"))
})


test_that(".noNAs works correctly", {
  expect_error(
    .noNAs(NA),
    "must not contain NA or NaN values"
  )
  expect_error(
    .noNAs(NaN, "argument"),
    "must not contain NA or NaN values"
  )
  expect_error(
    .noNAs(c(3, NA), "argument"),
    "must not contain NA or NaN values"
  )
  expect_error(
    .noNAs(c(3, NaN), "argument"),
    "must not contain NA or NaN values"
  )
  expect_null(.noNAs(3, "argument"))
  expect_null(.noNAs("foo", "argument"))
  expect_null(.noNAs(NULL, "argument"))
  expect_null(.noNAs(logical(0), "argument"))
  expect_null(.noNAs(numeric(0), "argument"))
  expect_null(.noNAs(character(0), "argument"))
  expect_null(.noNAs(Inf, "argument"))
  expect_null(.noNAs(TRUE, "argument"))
})

test_that(".MUST.isFinite works correctly", {

  expect_error(
    .MUST.isFinite(NULL, "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(NA, "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite("foo", "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(TRUE, "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(character(0), "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(numeric(0), "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(NaN, "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(Inf, "argument"),
    "must contain finite values"
  )
  expect_error(
    .MUST.isFinite(c(2, Inf), "argument"),
    "must contain finite values"
  )
  expect_null(.MUST.isFinite(3, "argument"))
  expect_null(.MUST.isFinite(1:4, "argument"))
})


# Type Checks -------------------------------------------------------------



test_that(".MUST.isLogical works correctly", {
  expect_error(
    .MUST.isLogical(NaN, "argument"),
    paste("must be of type", dQuote("logical"))
  )
  expect_error(
    .MUST.isLogical(3, "argument"),
    paste("must be of type", dQuote("logical"))
  )
  expect_error(
    .MUST.isLogical(1, "argument"),
    paste("must be of type", dQuote("logical"))
  )
  expect_error(
    .MUST.isLogical(0, "argument"),
    paste("must be of type", dQuote("logical"))
  )
  expect_error(
    .MUST.isLogical("TRUE", "argument"),
    paste("must be of type", dQuote("logical"))
  )
  expect_error(
    .MUST.isLogical(NULL, "argument"),
    paste("must be of type", dQuote("logical"))
  )
  expect_error(
    .MUST.isLogical(NA, "argument"),
    "must not contain NA or NaN values"
  )
  expect_null(.MUST.isLogical(TRUE, "argument"))
  expect_null(.MUST.isLogical(FALSE, "argument"))
  expect_null(.MUST.isLogical(logical(0), "argument"))
  expect_null(.MUST.isLogical(c(TRUE, FALSE), "argument"))
})

test_that(".MUST.isChar works correctly", {
  expect_error(
    .MUST.isChar(TRUE, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_error(
    .MUST.isChar(3, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_error(
    .MUST.isChar(1, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_error(
    .MUST.isChar(0, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_error(
    .MUST.isChar(NaN, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_error(
    .MUST.isChar(NA, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_error(
    .MUST.isChar(NULL, "argument"),
    paste("must be of type", dQuote("character"))
  )
  expect_null(.MUST.isChar("foo", "argument"))
  expect_null(.MUST.isChar(character(0), "argument"))
  expect_null(.MUST.isChar(c("foo", "bar"), "argument"))
})

test_that(".MUST.isNumeric works correctly", {

  expect_error(
    .MUST.isNumeric(NULL, "argument"),
    paste("must be of type", dQuote("numeric"))
  )
  expect_error(
    .MUST.isNumeric(NA, "argument"),
    paste("must be of type", dQuote("numeric"))
  )
  expect_error(
    .MUST.isNumeric("foo", "argument"),
    paste("must be of type", dQuote("numeric"))
  )
  expect_error(
    .MUST.isNumeric(TRUE, "argument"),
    paste("must be of type", dQuote("numeric"))
  )
  expect_null(.MUST.isNumeric(3, "argument"))
  expect_null(.MUST.isNumeric(numeric(0), "argument"))
  expect_null(.MUST.isNumeric(NaN, "argument"))
  expect_null(.MUST.isNumeric(Inf, "argument"))
  expect_null(.MUST.isNumeric(1:4, "argument"))
})

test_that(".MUST.isCharOrNumeric works correctly", {
  expect_error(
    .MUST.isCharOrNumeric(TRUE, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("numeric"))
  )
  expect_error(
    .MUST.isCharOrNumeric(NULL, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("numeric"))
  )
  expect_error(
    .MUST.isCharOrNumeric(NA, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("numeric"))
  )
  expect_error(
    .MUST.isCharOrNumeric(list(c("fee", "fie"), 1:2), "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("numeric"))
  )
  expect_null(.MUST.isCharOrNumeric("foo", "argument"))
  expect_null(.MUST.isCharOrNumeric(character(0), "argument"))
  expect_null(.MUST.isCharOrNumeric(3, "argument"))
  expect_null(.MUST.isCharOrNumeric(numeric(0), "argument"))
  expect_null(.MUST.isCharOrNumeric(NaN, "argument"))
  expect_null(.MUST.isCharOrNumeric(Inf, "argument"))
  expect_null(.MUST.isCharOrNumeric(1:4, "argument"))
  expect_null(.MUST.isCharOrNumeric(c("foo", "bar"), "argument"))
})

test_that(".MUST.isCharOrLogical works correctly", {
  expect_error(
    .MUST.isCharOrLogical(NULL, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(list(c("fee", "fie"), 1:2), "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(3, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(1, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(0, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(NaN, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(Inf, "argument"),
    paste("must be of type", dQuote("character"), "or", dQuote("logical"))
  )
  expect_error(
    .MUST.isCharOrLogical(NA, "argument"),
    "must not contain NA or NaN values"
  )
  expect_null(.MUST.isCharOrLogical("foo", "argument"))
  expect_null(.MUST.isCharOrLogical(character(0), "argument"))
  expect_null(.MUST.isCharOrLogical(c("foo", "bar"), "argument"))
  expect_null(.MUST.isCharOrLogical(TRUE, "argument"))
  expect_null(.MUST.isCharOrLogical(FALSE, "argument"))
  expect_null(.MUST.isCharOrLogical(logical(0), "argument"))
  expect_null(.MUST.isCharOrLogical(c(TRUE, FALSE), "argument"))
})


# Scalar Types ------------------------------------------------------------



test_that(".MUST.isTF works correctly", {
  expect_error(
    .MUST.isTF(NA, "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(NaN, "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(3, "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(1, "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(0, "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF("TRUE", "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(logical(0), "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(NULL, "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_error(
    .MUST.isTF(c(TRUE, FALSE), "argument"),
    paste("must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE"))
  )
  expect_null(.MUST.isTF(TRUE, "argument"))
  expect_null(.MUST.isTF(FALSE, "argument"))
})

test_that(".checkTF works", {
  expect_true(.checkTF(TRUE, default = TRUE))
  expect_false(.checkTF(FALSE, default = FALSE))
  expect_warning(.checkTF(as.logical(NA), default = TRUE))
  expect_warning(.checkTF(logical(0), default = TRUE))
  expect_warning(.checkTF(c(TRUE, FALSE), default = TRUE))
  expect_warning(.checkTF(NULL, default = TRUE))
  expect_warning(.checkTF(NA, default = TRUE))
  expect_warning(.checkTF(3, default = TRUE))
  expect_warning(.checkTF(1:3, default = TRUE))
})

test_that(".MUST.isTFOrAuto works correctly", {
  expect_error(
    .MUST.isTFOrAuto(NA, "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(NaN, "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(3, "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(1, "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(0, "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto("TRUE", "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto("AUTO", "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(logical(0), "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(character(0), "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(numeric(0), "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(c(TRUE, FALSE), "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_error(
    .MUST.isTFOrAuto(NULL, "argument"),
    paste("must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto"))
  )
  expect_null(.MUST.isTFOrAuto(TRUE, "argument"))
  expect_null(.MUST.isTFOrAuto(FALSE, "argument"))
  expect_null(.MUST.isTFOrAuto("auto", "argument"))
})

test_that(".MUST.isNumericScalar works", {
  expect_null(.MUST.isNumericScalar(3))
  expect_null(.MUST.isNumericScalar(pi))
  expect_null(.MUST.isNumericScalar(1.0))
  expect_error(.MUST.isNumericScalar("foo"), "must be a finite numeric scalar")
  expect_error(.MUST.isNumericScalar(NaN), "must be a finite numeric scalar")
  expect_error(.MUST.isNumericScalar(Inf), "must be a finite numeric scalar")
  expect_error(.MUST.isNumericScalar(NA), "must be a finite numeric scalar")
  expect_error(.MUST.isNumericScalar(numeric(0)),
               "must be a finite numeric scalar"
  )
  expect_error(.MUST.isNumericScalar(1:2), "must be a finite numeric scalar")
})

test_that(".MUST.isNonneg works correctly", {

  expect_error(
    .MUST.isNonneg(NULL, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(NA, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg("foo", "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(TRUE, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(character(0), "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(numeric(0), "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(1:2, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(NaN, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(Inf, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(-1, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_error(
    .MUST.isNonneg(-1.5, "argument"),
    "must be a finite and nonnegative scalar"
  )
  expect_null(.MUST.isNonneg(3, "argument"))
  expect_null(.MUST.isNonneg(0, "argument"))
})

test_that(".MUST.isPos works correctly", {

  expect_error(
    .MUST.isPos(NULL, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(NA, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos("foo", "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(TRUE, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(character(0), "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(numeric(0), "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(1:2, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(NaN, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(Inf, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(-1, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(0, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_error(
    .MUST.isPos(-1.5, "argument"),
    "must be a finite and strictly positive scalar"
  )
  expect_null(.MUST.isPos(3, "argument"))
})

test_that(".MUST.isInt works correctly", {

  expect_error(
    .MUST.isInt(NULL, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(NA, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt("foo", "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(TRUE, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(character(0), "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(numeric(0), "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(1:2, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(NaN, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(Inf, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(-1.5, "argument"),
    "must be a finite integer"
  )
  expect_error(
    .MUST.isInt(pi, "argument"),
    "must be a finite integer"
  )
  expect_null(.MUST.isInt(3, "argument"))
  expect_null(.MUST.isInt(3.0, "argument"))
})


test_that(".MUST.isNonnegInt works correctly", {

  expect_error(
    .MUST.isNonnegInt(NULL, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(NA, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt("foo", "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(TRUE, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(character(0), "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(numeric(0), "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(1:2, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(NaN, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(Inf, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(1.5, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(-1.5, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_error(
    .MUST.isNonnegInt(-1, "argument"),
    "must be a finite, nonnegative integer"
  )
  expect_null(.MUST.isNonnegInt(0, "argument"))
  expect_null(.MUST.isNonnegInt(3, "argument"))
  expect_null(.MUST.isNonnegInt(3.0, "argument"))
})

test_that(".MUST.isPosInt works correctly", {

  expect_error(
    .MUST.isPosInt(NULL, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(NA, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt("foo", "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(TRUE, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(character(0), "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(numeric(0), "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(1:2, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(NaN, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(Inf, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(1.5, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(-1.5, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(-1, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isPosInt(0, "argument"),
    "must be a finite, strictly positive integer"
  )
  expect_null(.MUST.isPosInt(3, "argument"))
  expect_null(.MUST.isPosInt(3.0, "argument"))
})

test_that(".MUST.isString works correctly", {
  expect_error(
    .MUST.isString(NULL, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(TRUE, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(3, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(1, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(0, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(NaN, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(NA, "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(c("foo", "bar"), "argument"),
    "must be a character string"
  )
  expect_error(
    .MUST.isString(character(0), "argument"),
    "must be a character string"
  )
  expect_null(.MUST.isString("foo", "argument"))
  expect_null(.MUST.isString("", "argument"))
})

test_that(".MUST.isNonemptyString works correctly", {
  expect_error(
    .MUST.isNonemptyString(NULL, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(TRUE, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(3, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(1, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(0, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(NaN, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(NA, "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(c("foo", "bar"), "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString(character(0), "argument"),
    "must be a nonempty character string"
  )
  expect_error(
    .MUST.isNonemptyString("", "argument"),
    "must be a nonempty character string"
  )
  expect_null(.MUST.isNonemptyString("foo", "argument"))
})


test_that(".MUST.isCharOrNumericScalar works correctly", {
  expect_error(
    .MUST.isCharOrNumericScalar(NA, "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(TRUE, "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(NULL, "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(list(c("fee", "fie"), 1:2), "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(1:2, "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(c("fee", "fie"), "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(character(0), "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(numeric(0), "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(NaN, "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_error(
    .MUST.isCharOrNumericScalar(Inf, "argument"),
    "must be a character string or finite numeric scalar"
  )
  expect_null(.MUST.isCharOrNumericScalar("foo", "argument"))
  expect_null(.MUST.isCharOrNumericScalar(3, "argument"))
  expect_null(.MUST.isCharOrNumericScalar(0, "argument"))
  expect_null(.MUST.isCharOrNumericScalar(pi, "argument"))
  expect_null(.MUST.isCharOrNumericScalar(-pi, "argument"))
})


test_that(".MUST.isStringOrInt works correctly", {
  expect_error(
    .MUST.isStringOrInt(NA, "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(TRUE, "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(NULL, "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(list(c("fee", "fie"), 1:2), "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(1:2, "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(c("fee", "fie"), "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(character(0), "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(numeric(0), "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(NaN, "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(Inf, "argument"),
    "must be a character string or finite integer"
  )
  expect_error(
    .MUST.isStringOrInt(pi, "argument"),
    "must be a character string or finite integer"
  )
  expect_null(.MUST.isStringOrInt("foo", "argument"))
  expect_null(.MUST.isStringOrInt(3, "argument"))
  expect_null(.MUST.isStringOrInt(-3, "argument"))
  expect_null(.MUST.isStringOrInt(0, "argument"))
})

test_that(".MUST.isStringOrPosInt works correctly", {
  expect_error(
    .MUST.isStringOrPosInt(NA, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(TRUE, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(NULL, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(list(c("fee", "fie"), 1:2), "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(1:2, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(c("fee", "fie"), "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(character(0), "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(numeric(0), "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(NaN, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(Inf, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(pi, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(-3, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_error(
    .MUST.isStringOrPosInt(0, "argument"),
    "must be a character string or a finite, strictly positive integer"
  )
  expect_null(.MUST.isStringOrPosInt("foo", "argument"))
  expect_null(.MUST.isStringOrPosInt(3, "argument"))
})


# Vector Types ------------------------------------------------------------

test_that(".isLogicalVector works", {
  expect_false(
    .isLogicalVector(NULL)
  )
  expect_true(
    .isLogicalVector(TRUE)
  )
  expect_false(
    .isLogicalVector(3)
  )
  expect_false(
    .isLogicalVector(NaN)
  )
  expect_true(
    .isLogicalVector(NA)
  )
  expect_false(
    .isLogicalVector(character(0))
  )
  expect_false(
    .isLogicalVector(numeric(0))
  )
  expect_false(
    .isLogicalVector(logical(0))
  )
  expect_false(
    .isLogicalVector(c("foo", NA))
  )
  expect_false(.isLogicalVector("foo"))
  expect_false(.isLogicalVector(c("foo", "bar")))
  expect_true(.isLogicalVector(c(TRUE, FALSE)))
})

test_that(".isNumericVector works", {
  expect_false(
    .isNumericVector(NULL)
  )
  expect_false(
    .isNumericVector(TRUE)
  )
  expect_true(
    .isNumericVector(3)
  )
  expect_true(
    .isNumericVector(pi)
  )
  expect_true(
    .isNumericVector(Inf)
  )
  expect_true(
    .isNumericVector(NaN)
  )
  expect_false(
    .isNumericVector(NA)
  )
  expect_true(
    .isNumericVector(c(NA, 3))
  )
  expect_true(
    .isNumericVector(c(Inf, 3))
  )
  expect_true(
    .isNumericVector(c(NaN, 3))
  )
  expect_false(
    .isNumericVector(character(0))
  )
  expect_false(
    .isNumericVector(numeric(0))
  )
  expect_false(
    .isNumericVector(logical(0))
  )
  expect_false(
    .isNumericVector(c("foo", NA))
  )
  expect_false(.isNumericVector("foo"))
  expect_false(.isNumericVector(c("foo", "bar")))
  expect_false(.isNumericVector(c(TRUE, FALSE)))
  expect_true(.isNumericVector(c(3, pi)))
  expect_true(.isNumericVector(1:4))
  expect_true(.isNumericVector(c(3, NA)))
  expect_false(.isNumericVector(as.factor(1:5)))
  expect_true(.isNumericVector(as.factor(1:5), factor_ok = TRUE))
  expect_true(.isNumericVector(as.factor(c("a", "b", "c")), factor_ok = TRUE))
})

test_that(".isIntegerVector works", {
  expect_false(
    .isIntegerVector(NULL)
  )
  expect_false(
    .isIntegerVector(TRUE)
  )
  expect_true(
    .isIntegerVector(3)
  )
  expect_false(
    .isIntegerVector(pi)
  )
  expect_false(
    .isIntegerVector(c(3, pi))
  )
  expect_false(
    .isIntegerVector(c(3, NaN))
  )
  expect_false(
    .isIntegerVector(c(3, Inf))
  )
  expect_false(
    .isIntegerVector(c(3, NA))
  )
  expect_false(
    .isIntegerVector(NaN)
  )
  expect_false(
    .isIntegerVector(NA)
  )
  expect_false(
    .isIntegerVector(character(0))
  )
  expect_false(
    .isIntegerVector(numeric(0))
  )
  expect_false(
    .isIntegerVector(logical(0))
  )
  expect_false(
    .isIntegerVector(c("foo", NA))
  )
  expect_false(.isIntegerVector("foo"))
  expect_false(.isIntegerVector(c("foo", "bar")))
  expect_false(.isIntegerVector(c(TRUE, FALSE)))
  expect_false(.isIntegerVector(c(3, pi)))
  expect_true(.isIntegerVector(c(-3.0, 2.0)))
  expect_true(.isIntegerVector(1:4))
  expect_false(.isIntegerVector(c(3, NA)))
  expect_false(.isIntegerVector(as.factor(1:5)))
  expect_true(.isIntegerVector(as.factor(1:5), factor_ok = TRUE))
  expect_true(.isIntegerVector(as.factor(c("a", "b", "c")), factor_ok = TRUE))
  expect_error(.MUST.isIntegerVector(c(3, NA)))
  expect_error(.MUST.isIntegerVector(as.factor(1:5)))
  expect_null(.MUST.isIntegerVector(as.factor(1:5), factor_ok = TRUE))
  expect_null(.MUST.isIntegerVector(as.factor(c("a", "b", "c")),
                                    factor_ok = TRUE)
  )
})

test_that(".isNonnegIntegerVector works", {
  expect_false(
    .isNonnegIntegerVector(NULL)
  )
  expect_false(
    .isNonnegIntegerVector(TRUE)
  )
  expect_true(
    .isNonnegIntegerVector(3)
  )
  expect_false(
    .isNonnegIntegerVector(pi)
  )
  expect_false(
    .isNonnegIntegerVector(-3)
  )
  expect_true(
    .isNonnegIntegerVector(0)
  )
  expect_false(
    .isNonnegIntegerVector(c(3, pi))
  )
  expect_false(
    .isNonnegIntegerVector(c(3, -3))
  )
  expect_false(
    .isNonnegIntegerVector(c(3, NaN))
  )
  expect_false(
    .isNonnegIntegerVector(c(3, Inf))
  )
  expect_false(
    .isNonnegIntegerVector(c(3, NA))
  )
  expect_false(
    .isNonnegIntegerVector(NaN)
  )
  expect_false(
    .isNonnegIntegerVector(NA)
  )
  expect_false(
    .isNonnegIntegerVector(character(0))
  )
  expect_false(
    .isNonnegIntegerVector(numeric(0))
  )
  expect_false(
    .isNonnegIntegerVector(logical(0))
  )
  expect_false(
    .isNonnegIntegerVector(c("foo", NA))
  )
  expect_false(.isNonnegIntegerVector("foo"))
  expect_false(.isNonnegIntegerVector(c("foo", "bar")))
  expect_false(.isNonnegIntegerVector(c(TRUE, FALSE)))
  expect_false(.isNonnegIntegerVector(c(3, pi)))
  expect_true(.isNonnegIntegerVector(1:4))
  expect_true(.isNonnegIntegerVector(0:4))
  expect_true(.isNonnegIntegerVector(c(0.0, 2.0)))
  expect_false(.isNonnegIntegerVector(-4:4))
  expect_false(.isNonnegIntegerVector(c(3, NA)))
  expect_false(.isNonnegIntegerVector(as.factor(1:5)))
  expect_true(.isNonnegIntegerVector(as.factor(1:5), factor_ok = TRUE))
  expect_true(.isNonnegIntegerVector(as.factor(c("a", "b", "c")),
                                     factor_ok = TRUE)
  )
})

test_that(".isPosIntegerVector works", {
  expect_false(
    .isPosIntegerVector(NULL)
  )
  expect_false(
    .isPosIntegerVector(TRUE)
  )
  expect_true(
    .isPosIntegerVector(3)
  )
  expect_false(
    .isPosIntegerVector(pi)
  )
  expect_false(
    .isPosIntegerVector(-3)
  )
  expect_false(
    .isPosIntegerVector(0)
  )
  expect_false(
    .isPosIntegerVector(c(3, pi))
  )
  expect_false(
    .isPosIntegerVector(c(3, -3))
  )
  expect_false(
    .isPosIntegerVector(c(3, NaN))
  )
  expect_false(
    .isPosIntegerVector(c(3, Inf))
  )
  expect_false(
    .isPosIntegerVector(c(3, NA))
  )
  expect_false(
    .isPosIntegerVector(NaN)
  )
  expect_false(
    .isPosIntegerVector(NA)
  )
  expect_false(
    .isPosIntegerVector(character(0))
  )
  expect_false(
    .isPosIntegerVector(numeric(0))
  )
  expect_false(
    .isPosIntegerVector(logical(0))
  )
  expect_false(
    .isPosIntegerVector(c("foo", NA))
  )
  expect_false(.isPosIntegerVector("foo"))
  expect_false(.isPosIntegerVector(c("foo", "bar")))
  expect_false(.isPosIntegerVector(c(TRUE, FALSE)))
  expect_false(.isPosIntegerVector(c(3, pi)))
  expect_true(.isPosIntegerVector(1:4))
  expect_true(.isPosIntegerVector(c(1.0, 2.0)))
  expect_false(.isPosIntegerVector(0:4))
  expect_false(.isPosIntegerVector(-4:4))
  expect_false(.isPosIntegerVector(c(3, NA)))
  expect_false(.isPosIntegerVector(as.factor(1:5)))
  expect_true(.isPosIntegerVector(as.factor(1:5), factor_ok = TRUE))
  expect_true(.isPosIntegerVector(as.factor(c("a", "b", "c")),
                                  factor_ok = TRUE)
  )
})

test_that(".MUST.isCharVector works correctly", {
  expect_error(
    .MUST.isCharVector(NULL, "argument"),
    "must be a nonempty character vector"
  )
  expect_error(
    .MUST.isCharVector(TRUE, "argument"),
    "must be a nonempty character vector"
  )
  expect_error(
    .MUST.isCharVector(3, "argument"),
    "must be a nonempty character vector"
  )
  expect_error(
    .MUST.isCharVector(NaN, "argument"),
    "must be a nonempty character vector"
  )
  expect_error(
    .MUST.isCharVector(NA, "argument"),
    "must be a nonempty character vector"
  )
  expect_error(
    .MUST.isCharVector(character(0), "argument"),
    "must be a nonempty character vector"
  )
  expect_null(
    .MUST.isCharVector(c("foo", NA), "argument")
  )
  expect_null(.MUST.isCharVector("foo", "argument"))
  expect_null(.MUST.isCharVector(c("foo", "bar"), "argument"))
  expect_error(.MUST.isCharVector(as.factor(c("foo", "bar"))))
  expect_error(.MUST.isCharVector(as.factor(NULL), factor_ok = TRUE))
  expect_null(.MUST.isCharVector(as.factor(c("foo", "bar")), factor_ok = TRUE))
})

test_that(".MUST.isCharOrNumericVector works correctly", {
  expect_error(
    .MUST.isCharOrNumericVector(NULL, "argument"),
    paste("must be a nonempty vector of type",
          dQuote("character"), "or", dQuote("numeric")
    )
  )
  expect_error(
    .MUST.isCharOrNumericVector(TRUE, "argument"),
    paste("must be a nonempty vector of type",
          dQuote("character"), "or", dQuote("numeric")
    )
  )
  expect_error(
    .MUST.isCharOrNumericVector(NA, "argument"),
    paste("must be a nonempty vector of type",
          dQuote("character"), "or", dQuote("numeric")
    )
  )
  expect_error(
    .MUST.isCharOrNumericVector(list(c("fee", "fie"), 1:2), "argument"),
    paste("must be a nonempty vector of type",
          dQuote("character"), "or", dQuote("numeric")
    )
  )
  expect_error(
    .MUST.isCharOrNumericVector(character(0), "argument"),
    paste("must be a nonempty vector of type",
          dQuote("character"), "or", dQuote("numeric")
    )
  )
  expect_error(
    .MUST.isCharOrNumericVector(numeric(0), "argument"),
    paste("must be a nonempty vector of type",
          dQuote("character"), "or", dQuote("numeric")
    )
  )
  expect_null(
    .MUST.isCharOrNumericVector(c(3, NA), "argument")
  )

  expect_null(
    .MUST.isCharOrNumericVector(NaN, "argument")
  )
  expect_null(
    .MUST.isCharOrNumericVector(c(2.1, NaN), "argument")
  )
  expect_null(
    .MUST.isCharOrNumericVector(c(2.1, NA), "argument")
  )
  expect_null(
    .MUST.isCharOrNumericVector(Inf, "argument")
  )
  expect_null(
    .MUST.isCharOrNumericVector(c(2.1, Inf), "argument")
  )
  expect_null(
    .MUST.isCharOrNumericVector(c("foo", NA), "argument")
  )
  expect_null(.MUST.isCharOrNumericVector("foo", "argument"))
  expect_null(.MUST.isCharOrNumericVector(c("foo", "bar"), "argument"))
  expect_null(.MUST.isCharOrNumericVector(3, "argument"))
  expect_null(.MUST.isCharOrNumericVector(1:4, "argument"))
  expect_error(.MUST.isCharOrNumericVector(as.factor(c("foo", "bar"))))
  expect_error(.MUST.isCharOrNumericVector(as.factor(NULL), factor_ok = TRUE))
  expect_null(.MUST.isCharOrNumericVector(as.factor(c("foo", "bar")),
                                          factor_ok = TRUE)
  )
})

test_that(".MUST.isCharOrIntegerVector works correctly", {
  expect_error(
    .MUST.isCharOrIntegerVector(NULL, "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )
  expect_error(
    .MUST.isCharOrIntegerVector(TRUE, "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )
  expect_error(
    .MUST.isCharOrIntegerVector(NA, "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )
  expect_error(
    .MUST.isCharOrIntegerVector(list(c("fee", "fie"), 1:2), "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )
  expect_error(
    .MUST.isCharOrIntegerVector(character(0), "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )
  expect_error(
    .MUST.isCharOrIntegerVector(numeric(0), "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )
  expect_error(
    .MUST.isCharOrIntegerVector(c(3, NA), "argument"),
    "must be a nonempty character vector or integer-valued vector"
  )

  expect_error(
    .MUST.isCharOrIntegerVector(NaN, "argument")
  )
  expect_error(
    .MUST.isCharOrIntegerVector(c(2, NaN), "argument")
  )
  expect_error(
    .MUST.isCharOrIntegerVector(c(2, NA), "argument")
  )
  expect_error(
    .MUST.isCharOrIntegerVector(Inf, "argument")
  )
  expect_error(
    .MUST.isCharOrIntegerVector(c(2, Inf), "argument")
  )
  expect_error(
    .MUST.isCharOrIntegerVector(c(2, pi), "argument")
  )
  expect_null(
    .MUST.isCharOrIntegerVector(c(-1.0, 3.000), "argument")
  )
  expect_null(
    .MUST.isCharOrIntegerVector(c("foo", NA), "argument")
  )
  expect_null(.MUST.isCharOrIntegerVector("foo", "argument"))
  expect_null(.MUST.isCharOrIntegerVector(c("foo", "bar"), "argument"))
  expect_null(.MUST.isCharOrIntegerVector(3, "argument"))
  expect_null(.MUST.isCharOrIntegerVector(1:4, "argument"))
  expect_error(.MUST.isCharOrIntegerVector(as.factor(c("foo", "bar"))))
  expect_error(.MUST.isCharOrIntegerVector(as.factor(NULL), factor_ok = TRUE))
  expect_null(.MUST.isCharOrIntegerVector(as.factor(c("foo", "bar")),
                                          factor_ok = TRUE)
  )
})


# Regular Expressions -----------------------------------------------------

test_that(".isValidFilenamePart works", {
  expect_true(.isValidFilenamePart("foo"))
  expect_true(.isValidFilenamePart("foo123"))
  expect_true(.isValidFilenamePart("foo_bar"))
  expect_true(.isValidFilenamePart("foo-bar"))
  expect_false(.isValidFilenamePart("foo bar"))
  expect_false(.isValidFilenamePart("-foo"))
  expect_false(.isValidFilenamePart("foo_"))
  expect_false(.isValidFilenamePart("_foo_"))
  expect_false(.isValidFilenamePart("foo!"))
  expect_false(.isValidFilenamePart(""))
})

test_that(".isValidFilename works", {

  expect_true(.isValidFilename("foo.txt"))
  expect_true(.isValidFilename("foo_bar.txt"))
  expect_true(.isValidFilename("foo-bar.txt"))

  expect_false(.isValidFilename("foo bar.txt"))
  expect_false(.isValidFilename("foo.bar.txt"))
  expect_false(.isValidFilename("_foo.txt"))
  expect_false(.isValidFilename("foo.txt_"))
  expect_false(.isValidFilename("foo."))
  expect_false(.isValidFilename(".foo"))
  expect_false(.isValidFilename(""))
})

test_that(".sanitizeFilenamePart works", {
  expect_equal(.sanitizeFilenamePart("-foo"), "foo")
  expect_equal(.sanitizeFilenamePart("foo"), "foo")
  expect_equal(.sanitizeFilenamePart("foo_"), "foo")
  expect_equal(.sanitizeFilenamePart("_foo_"), "foo")
  expect_equal(.sanitizeFilenamePart("foo!"), "foo")
  expect_equal(.sanitizeFilenamePart(".foo"), "foo")
  expect_equal(.sanitizeFilenamePart("foo-bar"), "foo-bar")
  expect_equal(.sanitizeFilenamePart("foo_bar"), "foo_bar")
  expect_equal(.sanitizeFilenamePart("foo bar"), "foo_bar")
  expect_equal(.sanitizeFilenamePart("foo/bar"), "foo_bar")
  expect_equal(.sanitizeFilenamePart("foo.bar"), "foo_bar")
  expect_equal(
    .sanitizeFilenamePart("--_____-_--__-_-f-_-o__---_-_o-___--__-___-_-__!"),
    "f-o__-_o"
  )
  expect_equal(
    .sanitizeFilenamePart("--_- -_ -  __--  -__-_foo-bar_baz--!bing.txt__"),
    "foo-bar_baz-_bing_txt"
  )
  expect_equal(.sanitizeFilenamePart("!-_."), "")
})

# Lists and Data Frames ---------------------------------------------------

test_that(".MUST.isNamedList works", {
  expect_null(.MUST.isNamedList(list(foo = 1:2, bar = diag(2))))
  expect_error(.MUST.isNamedList(list(1:2, diag(2))), "must be a named list")
})

test_that(".MUST.isDataFrame works", {
  expect_null(.MUST.isDataFrame(as.data.frame(diag(2))))
  expect_null(.MUST.isDataFrame(data.frame(foo = 1:2, bar = c("bar", "baz"))))
  expect_error(.MUST.isDataFrame(diag(2)),
               "could not be coerced to a data frame"
  )
})

test_that(".MUST.hasMultipleRows works", {
  expect_null(.MUST.hasMultipleRows(diag(2)))
  expect_null(.MUST.hasMultipleRows(as.data.frame(diag(2))))
  expect_null(
    .MUST.hasMultipleRows(data.frame(foo = 1:2, bar = c("bar", "baz")))
  )
  expect_error(.MUST.hasMultipleRows(diag(1)), "must contain at least two rows")
  expect_error(.MUST.hasMultipleRows(data.frame(foo = 1, bar = "bar")),
               "must contain at least two rows"
  )
  expect_error(.MUST.hasMultipleRows(matrix()), "must contain at least two rows")
  expect_error(
    .MUST.hasMultipleRows(data.frame(foo = numeric(0), bar = character(0))),
    "must contain at least two rows"
  )
  expect_error(.MUST.hasMultipleRows(NULL), "must contain at least two rows")
  foo <- data.frame(a = 1, b = 1, c = 2)
  foo2 <- rbind(foo, foo)
  expect_error(
    .MUST.hasMultipleRows(foo),
    "must contain at least two rows"
  )
  expect_null(.MUST.hasMultipleRows(foo2))
})



# Column References -------------------------------------------------------


test_that(".MUST.isDataColref works", {
  expect_null(.MUST.isDataColref("CloneSeq", dat, "argument"))
  expect_null(.MUST.isDataColref(1, dat, "argument"))
  expect_null(.orNull(.MUST.isDataColref, "CloneSeq", data = dat))
  expect_true(.orNull(.MUST.isDataColref, NULL, data = dat))
  expect_null(.MUST.isDataColrefs("CloneSeq", dat, "argument"))
  expect_null(
    .MUST.isDataColrefs(c("CloneSeq", "CloneCount"), dat, "argument")
  )
  expect_null(
    .orNull(.MUST.isDataColrefs, c("CloneSeq", "CloneCount"), data = dat)
  )
  expect_true(.orNull(.MUST.isDataColrefs, NULL, data = dat))
  expect_error(.MUST.isDataColref(NULL, dat, "argument"),
               "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(TRUE, dat, "argument"),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(NA, dat, "argument"),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(c(1, 2), dat, "argument"),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(c("CloneSeq", "CloneCount"), dat, "argument"),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(character(0), dat, "argument"),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(numeric(0), dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref("foo", dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(NaN, dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(Inf, dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(1.1, dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(10, dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(0, dat),
    "must specify a data column"
  )
  expect_error(
    .MUST.isDataColref(-5, dat),
    "must specify a data column"
  )

  expect_error(
    .MUST.isDataColrefs(NULL, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(TRUE, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(NA, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(character(0), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(numeric(0), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs("foo", dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(NaN, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(Inf, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(1.1, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(10, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(0, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(-5, dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c("foo", NA), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, NA), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, NaN), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c("CloneSeq", "foo"), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, Inf), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, 1.1), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, 10), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, 0), dat),
    "must specify one or more data columns"
  )
  expect_error(
    .MUST.isDataColrefs(c(1, -10), dat),
    "must specify one or more data columns"
  )

})


test_that(".checkDataColrefs works", {
  expect_equal(.checkDataColrefs("CloneSeq", dat), "CloneSeq")
  expect_equal(.checkDataColrefs(c("CloneSeq", "CloneCount"), dat),
               c("CloneSeq", "CloneCount")
  )
  expect_equal(.checkDataColrefs(colnames(dat), dat),
               colnames(dat)
  )
  expect_equal(.checkDataColrefs(1, dat), 1)
  expect_equal(.checkDataColrefs(1:4, dat), 1:4)
  expect_warning(.checkDataColrefs("foo", dat), "one or more values")
  expect_warning(.checkDataColrefs(c("CloneSeq", "foo"), dat),
                 "one or more values"
  )
  expect_warning(.checkDataColrefs(0, dat), "one or more values")
  expect_warning(.checkDataColrefs(-1, dat), "one or more values")
  expect_warning(.checkDataColrefs(20, dat), "one or more values")
  expect_warning(.checkDataColrefs(c(1, 20), dat), "one or more values")
  expect_warning(.checkDataColrefs(c(1, -1), dat), "one or more values")
  expect_warning(.checkDataColrefs(c(1, 0), dat), "one or more values")
  expect_warning(.checkDataColrefs(character(0), dat), "one or more values")
  expect_warning(.checkDataColrefs(numeric(0), dat), "one or more values")
  expect_warning(.checkDataColrefs(TRUE, dat), "one or more values")
})


test_that(".checkDataForColref works", {
  expect_equal(.checkDataForColref("CloneSeq", dat, .isCharVector), "CloneSeq")
  expect_equal(.checkDataForColref("CloneCount", dat, .isIntegerVector),
               "CloneCount"
  )
  expect_warning(.checkDataForColref("CloneSeq", dat, .isIntegerVector),
                 "with the correct properties."
  )
  expect_warning(.checkDataForColref("CloneCount", dat, .isCharVector),
                 "with the correct properties."
  )
})

# AIRR-Seq Data -----------------------------------------------------------



test_that(".isValidSeqVector works correctly", {
  expect_error(
    .MUST.isValidSeqVector(NULL),
    "must be coercible to a character vector with at least one non-NA value"
  )
  expect_error(
    .MUST.isValidSeqVector(character(0)),
    "must be coercible to a character vector with at least one non-NA value"
  )
  expect_error(
    .MUST.isValidSeqVector(NA),
    "must be coercible to a character vector with at least one non-NA value"
  )
  expect_error(
    .MUST.isValidSeqVector(c(NA, NA)),
    "must be coercible to a character vector with at least one non-NA value"
  )
  expect_null(.MUST.isValidSeqVector("foo"))
  expect_null(.MUST.isValidSeqVector(1))
  expect_null(.MUST.isValidSeqVector(c("foo", "bar")))
  expect_null(.MUST.isValidSeqVector(c("foo", NA)))
  expect_null(.MUST.isValidSeqVector(1:2))
  expect_null(.MUST.isValidSeqVector(c(1, NA, NaN, Inf)))
  expect_null(.MUST.isValidSeqVector(as.factor(c("foo", "bar", "bar"))))
})

test_that(".MUST.isSeqColref works correctly", {
  expect_error(
    .MUST.isSeqColref(NULL, dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColrefs(c("CloneSeq", "CloneCount", "CloneFreq"), dat),
    "does not reference one or two columns of"
  )
  expect_error(
    .MUST.isSeqColref(character(0), dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref(numeric(0), dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref("foo", dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref(1.1, dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref(0, dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref(10, dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref(-5, dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColref(c("CloneSeq", "CloneCount"), dat),
    "does not reference a column of"
  )
  expect_error(
    .MUST.isSeqColrefs(c("CloneSeq", "foo"), dat),
    "does not reference one or two columns of"
  )
  expect_error(
    .MUST.isSeqColrefs(c(-1, -5), dat),
    "does not reference one or two columns of"
  )
  expect_error(
    .MUST.isSeqColrefs(c(1, -5), dat),
    "does not reference one or two columns of"
  )
  expect_error(
    .MUST.isSeqColrefs(c(1, 1.1), dat),
    "does not reference one or two columns of"
  )
  expect_null(.MUST.isSeqColref("CloneSeq", dat))
  expect_null(.MUST.isSeqColref("CloneCount", dat))
  expect_null(.MUST.isSeqColrefs(c("CloneSeq", "CloneCount"), dat))
})


test_that(".isCountColref works", {
  expect_null(.MUST.isCountColref("CloneCount", dat))
  expect_null(.MUST.isCountColref("CloneFrequency", dat))
  expect_error(.MUST.isCountColref("CloneSeq", dat),
               "that contains numeric values"
  )
})

test_that(".checkCountCol works", {
  expect_equal(.checkCountCol("CloneCount", dat), "CloneCount")
  expect_equal(.checkCountCol("CloneFrequency", dat), "CloneFrequency")
  expect_warning(.checkCountCol("CloneSeq", dat),
                 "or does not specify a numeric variable of"
  )
})

# Network Objects ---------------------------------------------------------


test_that("network object type checks work correctly", {

  expect_null(.MUST.isIgraph(net$igraph))
  expect_null(.MUST.isIgraph(net0$igraph))
  expect_error(
    .MUST.isIgraph(net$node_data, "argument"),
    paste("must be of class", dQuote("igraph"))
  )

  expect_null(.MUST.isGgraph(net$plots[[1]]))
  expect_null(.MUST.isGgraph(net0$plots[[1]]))
  expect_error(
    .MUST.isGgraph(net$plots, "argument"),
    paste("must be of class", dQuote("ggraph"))
  )

  expect_null(
    .MUST.isAdjacencyMatrix(matrix(0, nrow = 2, ncol = 2), "argument")
  )
  expect_error(
    .MUST.isAdjacencyMatrix(net, "argument"),
    "must be a symmetric matrix"
  )
  expect_error(
    .MUST.isAdjacencyMatrix(net, "argument"),
    "must be a symmetric matrix"
  )
  expect_error(
    .MUST.isAdjacencyMatrix(matrix(0, nrow = 1, ncol = 2), "argument"),
    "must be a symmetric matrix"
  )
  expect_error(
    .MUST.isAdjacencyMatrix(matrix(2, nrow = 2, ncol = 2), "argument"),
    "must be a symmetric matrix"
  )

  expect_null(.MUST.isLayout(net$plots$graph_layout))
  expect_null(.MUST.isLayout(net0$plots$graph_layout))
  expect_error(.MUST.isLayout(diag(3)),
               "must be a nonempty two-column numeric matrix"
  )

  expect_null(.MUST.isPlotlist(net$plots))
  expect_null(.MUST.isPlotlist(net0$plots))
  expect_error(.MUST.isPlotlist(net), "and possibly a two-column numeric")

})


# Correspondence of network objects ---------------------------------------



test_that("Checks for correspondence of network objects work", {

  expect_null(.MUST.doesIgraphMatchData(net$igraph, net$node_data))
  expect_null(.MUST.doesIgraphMatchData(net0$igraph, net0$node_data))
  expect_error(
    .MUST.doesIgraphMatchData(net$igraph, net0$node_data),
    "must have length equal to row dimension of"
  )

  expect_null(.MUST.doesIgraphMatchMatrix(net$igraph, net$adjacency_matrix))
  expect_null(.MUST.doesIgraphMatchMatrix(net0$igraph, net0$adjacency_matrix))
  expect_error(
    .MUST.doesIgraphMatchMatrix(net$igraph, net0$adjacency_matrix),
    "must have length equal to row dimension of"
  )

  expect_null(.MUST.doesDataMatchMatrix(net$node_data, net$adjacency_matrix))
  expect_null(.MUST.doesDataMatchMatrix(net0$node_data, net0$adjacency_matrix))
  expect_error(
    .MUST.doesDataMatchMatrix(net$node_data, net0$adjacency_matrix),
    "must have the same row dimension"
  )

  expect_null(.MUST.doesPlotMatchData(net$plots[[1]], net$node_data))
  expect_null(.MUST.doesPlotMatchData(net0$plots[[1]], net0$node_data))
  expect_error(
    .MUST.doesPlotMatchData(net$plots[[1]], net0$node_data),
    "has node count not equal to the row dimension of"
  )

  expect_null(.MUST.doesLayoutMatchData(net$plots$graph_layout, net$node_data))
  expect_null(
    .MUST.doesLayoutMatchData(net0$plots$graph_layout, net0$node_data)
  )
  expect_error(
    .MUST.doesLayoutMatchData(net$plots$graph_layout, net0$node_data),
    "must be NULL or a two-column matrix"
  )


  # Network list ------------------------------------------------------------


})

test_that("Checks for network output list work", {

  expect_null(.MUST.hasIgraph(net))
  expect_null(.MUST.hasIgraph(net0))
  expect_error(.MUST.hasIgraph(net$plots), "must contain an object of class")

  expect_null(.MUST.hasAdjacencyMatrix(net))
  expect_null(.MUST.hasAdjacencyMatrix(net0))
  expect_error(.MUST.hasAdjacencyMatrix(net$plots), "must contain a symmetric")

  expect_null(.MUST.hasNodeData(net, "argument"))
  expect_null(.MUST.hasNodeData(net0, "argument"))
  expect_error(.MUST.hasNodeData(net$plots), "nonempty data frame named")

  expect_null(.MUST.hasDetails(net))
  expect_null(.MUST.hasDetails(net0))
  expect_error(.MUST.hasDetails(net$plots), "nonempty named list named")

  expect_null(.MUST.hasClusterData(net, "argument"))
  expect_error(.MUST.hasClusterData(net0), "integer-valued variable named")

  expect_null(.MUST.hasPlots(net))
  expect_null(.MUST.hasPlots(net0))
  expect_error(.MUST.hasPlots(net$plots), "possibly a two-column")

  expect_null(.MUST.isBaseNetworkOutput(net, "argument"))
  expect_null(.MUST.isBaseNetworkOutput(net0, "argument"))
  expect_error(
    .MUST.isBaseNetworkOutput(net$node_data, "argument"),
    paste("must be a named list containing elements",
          paste0(dQuote("igraph"), ","),
          dQuote("adjacency_matrix"), "and", paste0(dQuote("node_data"), ","),
          "all corresponding to the same network"
    )
  )
})



# File Input Arguments ----------------------------------------------------


test_that(".isInputType works correctly", {
  valid_input_types <- c("csv", "table", "tsv", "txt", "rds", "rda")
  for (i in 1:length(valid_input_types)) {
    expect_true(.isInputType(valid_input_types[[i]]))
  }
  expect_error(
    .MUST.isInputType("foo"),
    "must be one of:"
  )

})


test_that(".checkargs.InputFiles works correctly", {
  x <- NULL
  y <- NULL
  z <- NULL
  file_a <- file.path(tempdir(), "a.rds")
  file_b <- file.path(tempdir(), "b.rds")
  file_c <- file.path(tempdir(), "c.rds")
  con_b <- file(file_b)
  con_c <- gzfile(file_c)
  saveRDS(x, file = file_a)
  saveRDS(y, file = file_b)
  saveRDS(z, file = file_c)
  good_file_list <- c(file_a, file_b, file_c)
  good_file_list_b <- list(file_a, con_b, con_c)
  bad_file_list <- c(file_a, file_b, "foo")
  bad_file_list_b <- list(file_a, con_b, "foo")
  bad_file_list_c <- c(file_a, file_a, file_b)
  bad_file_list_d <- list(file_a, con_b, diag(3))
  bad_file_list_e <- list(file_a, con_b, c("foo", "bar"))
  bad_file_list_f <- 5
  good_data_symbols_1 <- c("x", "y", "z")
  good_data_symbols_2 <- "x"
  good_readargs <- list(sep = "", header = TRUE, row.names = 1)
  bad_readargs <- c(sep = "", header = TRUE)
  bad_readargs_b <- list(sep = "", foo = "bar")
  expect_null(
    .checkargs.InputFiles(good_file_list, "rda", good_data_symbols_1, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "rda", good_data_symbols_2, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "rds", good_data_symbols_2, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "csv", good_data_symbols_2, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "rds", good_data_symbols_2, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "txt", good_data_symbols_2, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "tsv", good_data_symbols_2, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "tsv", NULL, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "table", NULL, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "table", 3, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list_b, "tsv", NULL, TRUE, "")
  )
  expect_null(
    .checkargs.InputFiles(good_file_list, "tsv", NULL, TRUE, "", good_readargs)
  )
  expect_error(
    .checkargs.InputFiles(good_file_list, "tsv", NULL, TRUE, "", bad_readargs),
    "must be a named list"
  )
  expect_error(
    .checkargs.InputFiles(good_file_list, "tsv", NULL, "", ""),
    "must evaluate to"
  )
  expect_error(
    .checkargs.InputFiles(good_file_list, "tsv", NULL, TRUE, TRUE),
    "must be a character string"
  )
  expect_error(
    .checkargs.InputFiles(good_file_list, "tsv", NULL, TRUE, "", bad_readargs),
    "must be a named list"
  )
  expect_warning(.checkReadArgs(bad_readargs_b, TRUE, ""), "dropping argument")
  expect_error(
    .checkargs.InputFiles(good_file_list, "rda", 3, TRUE, ""),
    "must be a nonempty character vector"
  )
  expect_error(
    .checkargs.InputFiles(good_file_list, "rda", c("x", "y"), TRUE, ""),
    "must have length 1 or equal to that of"
  )
  expect_error(
    .checkargs.InputFiles(bad_file_list, "rds", NULL, TRUE, ""),
    "nonexistent files"
  )
  expect_error(
    .checkargs.InputFiles(bad_file_list_b, "rds", NULL, TRUE, ""),
    "nonexistent files"
  )
  expect_error(
    .checkargs.InputFiles(bad_file_list_c, "rds", NULL, TRUE, ""),
    "duplicate values"
  )
  expect_error(
    .checkargs.InputFiles(bad_file_list_d, "rds", NULL, TRUE, ""),
    "contains elements other than"
  )
  expect_error(
    .checkargs.InputFiles(bad_file_list_e, "rds", NULL, TRUE, ""),
    "contains elements other than"
  )
  expect_error(
    .checkargs.InputFiles(bad_file_list_f, "rds", NULL, TRUE, ""),
    "or a list of character strings and connections"
  )
})


# File Output Arguments ---------------------------------------------------

test_that(".isOutputType works correctly", {
  expect_false(
    .isOutputType("csv")
  )
  expect_false(
    .isOutputType("tsv", "findPublicClusters")
  )
  expect_false(
    .isOutputType("individual", "findPublicClusters")
  )
  expect_false(
    .isOutputType("individual", "findAssociatedClones")
  )
  expect_true(
    .isOutputType("table", "findAssociatedClones")
  )
  expect_false(
    .isOutputType("individual", "generic")
  )
  expect_true(
    .isOutputType("individual")
  )
  expect_true(
    .isOutputType("rds", "findPublicClusters")
  )
  expect_true(
    .isOutputType("tsv", "findAssociatedClones")
  )
  expect_true(
    .isOutputType("table", "generic")
  )
})

test_that(".checkOutputName works", {
  expect_equal(.checkOutputName("foo_bar-123"), "foo_bar-123")
  expect_warning(.checkOutputName(NULL),
                 paste("Defaulting to", dQuote("MyRepSeqNetwork"))
  )
  expect_warning(.checkOutputName(c("foo", "bar")),
                 paste("Defaulting to", dQuote("MyRepSeqNetwork"))
  )
  expect_warning(.checkOutputName(""),
                 paste("Defaulting to", dQuote("MyRepSeqNetwork"))
  )
  expect_warning(.checkOutputName("_-."),
                 paste("Value changed to", dQuote("MyRepSeqNetwork"))
  )
  expect_warning(.checkOutputName("foo.bar"),
                 paste("Value changed to", dQuote("foo_bar"))
  )
  expect_warning(.checkOutputName("_foo-bar."),
                 paste("Value changed to", dQuote("foo-bar"))
  )
})

test_that(".checkOutfileLayout works", {
  expect_equal(.checkOutfileLayout("foo", net$plots), "foo")
  expect_equal(.checkOutfileLayout("foo", net0$plots), "foo")
  expect_warning(.checkOutfileLayout("foo", net1), "valid layout matrix")
})


# Network Analysis Arguments ----------------------------------------------



test_that(".isDistType works correctly", {
  expect_error(
    .MUST.isDistType("foo"),
    paste("must be", dQuote("hamming"), "or", dQuote("levenshtein"))
  )
  expect_null(.MUST.isDistType("levenshtein"))
  expect_null(.MUST.isDistType("lev"))
  expect_null(.MUST.isDistType("l"))
  expect_null(.MUST.isDistType("Levenshtein"))
  expect_null(.MUST.isDistType("Lev"))
  expect_null(.MUST.isDistType("L"))
  expect_null(.MUST.isDistType("h"))
  expect_null(.MUST.isDistType("ham"))
  expect_null(.MUST.isDistType("hamming"))
  expect_null(.MUST.isDistType("H"))
  expect_null(.MUST.isDistType("Ham"))
  expect_null(.MUST.isDistType("Hamming"))

  expect_equal(.matchDistType("levenshtein"), "levenshtein")
  expect_equal(.matchDistType("lev"), "levenshtein")
  expect_equal(.matchDistType("l"), "levenshtein")
  expect_equal(.matchDistType("Levenshtein"), "levenshtein")
  expect_equal(.matchDistType("Lev"), "levenshtein")
  expect_equal(.matchDistType("L"), "levenshtein")
  expect_equal(.matchDistType("h"), "hamming")
  expect_equal(.matchDistType("ham"), "hamming")
  expect_equal(.matchDistType("hamming"), "hamming")
  expect_equal(.matchDistType("H"), "hamming")
  expect_equal(.matchDistType("Ham"), "hamming")
  expect_equal(.matchDistType("Hamming"), "hamming")

  expect_equal(.checkDistType("levenshtein"), "levenshtein")
  expect_equal(.checkDistType("lev"), "levenshtein")
  expect_equal(.checkDistType("l"), "levenshtein")
  expect_equal(.checkDistType("Levenshtein"), "levenshtein")
  expect_equal(.checkDistType("Lev"), "levenshtein")
  expect_equal(.checkDistType("L"), "levenshtein")
  expect_equal(.checkDistType("h"), "hamming")
  expect_equal(.checkDistType("ham"), "hamming")
  expect_equal(.checkDistType("hamming"), "hamming")
  expect_equal(.checkDistType("H"), "hamming")
  expect_equal(.checkDistType("Ham"), "hamming")
  expect_equal(.checkDistType("Hamming"), "hamming")
  expect_warning(.checkDistType("foo"), "invalid.")
})


test_that(".checkStatsToInclude works correctly", {
  expect_warning(
    .checkStatsToInclude(NULL)
  )
  expect_warning(
    .checkStatsToInclude("foo")
  )
  expect_warning(
    .checkStatsToInclude(c("degree" = TRUE, "cluster_id" = TRUE))
  )
  expect_equal(
    .checkStatsToInclude(chooseNodeStats()), chooseNodeStats()
  )
  expect_equal(
    .checkStatsToInclude(exclusiveNodeStats()), exclusiveNodeStats()
  )
  expect_equal(.checkStatsToInclude("all"), chooseNodeStats(all_stats = TRUE))
})


test_that(".isClusterFun works correctly", {
  expect_false(
    .isClusterFun(NULL)
  )
  expect_false(
    .isClusterFun(diag(5))
  )
  expect_false(
    .isClusterFun(log)
  )
  expect_false(
    .isClusterFun(NA)
  )
  expect_false(
    .isClusterFun(0)
  )
  expect_false(
    .isClusterFun(c("cluster_fast_greedy", "cluster_leiden"))
  )
  expect_false(
    .isClusterFun("log")
  )
  expect_true(
    .isClusterFun("cluster_fast_greedy")
  )
  expect_true(
    .isClusterFun("cluster_leiden")
  )
  expect_true(
    .isClusterFun("fast_greedy")
  )
})


# Plotting Arguments ------------------------------------------------------


test_that(".checkSizeNodesBy works correctly", {
  dat <- simulateToyData(sample_size = 5)
  expect_warning(
    .checkSizeNodesBy(1:2, dat)
  )
  expect_warning(
    .checkSizeNodesBy("foo", dat)
  )
  expect_warning(
    .checkSizeNodesBy(0, dat)
  )
  expect_equal(
    .checkSizeNodesBy("CloneSeq", dat), "CloneSeq"
  )
  expect_equal(
    .checkSizeNodesBy("CloneCount", dat), "CloneCount"
  )
  expect_equal(
    .checkSizeNodesBy(1.5, dat), 1.5
  )
  expect_equal(
    .checkSizeNodesBy(NULL, dat), NULL
  )
})


test_that(".checkColorNodesBy works correctly", {
  dat <- simulateToyData(sample_size = 5)
  expect_warning(.checkColorNodesBy("foo", dat))
  expect_warning(.checkColorNodesBy("degree", dat))
  expect_warning(.checkColorNodesBy(c("CloneSeq", "degree"), dat))
  expect_equal(.checkColorNodesBy(c("CloneSeq", "CloneCount"), dat),
               c("CloneSeq", "CloneCount"),
               ignore_attr = TRUE
  )
  expect_equal(.checkColorNodesBy("CloneSeq", dat), "CloneSeq",
               ignore_attr = TRUE
  )
  expect_warning(.checkColorNodesBy("auto", dat))
  expect_equal(.checkColorNodesBy("auto", dat, auto_ok = TRUE),
               "auto", ignore_attr = TRUE)
  expect_equal(.checkColorNodesBy("degree", dat, node_stats = TRUE),
               "degree", ignore_attr = TRUE
  )
  expect_equal(
    .checkColorNodesBy(c("CloneSeq", "degree"), dat, node_stats = TRUE),
    c("CloneSeq", "degree"), ignore_attr = TRUE
  )
  expect_equal(.checkColorNodesBy("foo", dat, plots = FALSE), "foo",
               ignore_attr = TRUE
  )
  expect_equal(.checkColorNodesBy(NULL, dat), NULL, ignore_attr = TRUE)
})

test_that(".checkColorScheme works correctly", {
  expect_warning(.checkColorScheme("foo", "CloneSeq")
  )
  expect_warning(
    .checkColorScheme(c("viridis", "foo"), c("CloneSeq", "CloneCount"))
  )
  expect_warning(
    .checkColorScheme(c("foo", "bar", "baz"), c("CloneSeq", "CloneCount")),
  )
  expect_equal(
    .checkColorScheme("viridis", c("CloneSeq", "CloneCount")),
    "viridis", ignore_attr = TRUE
  )
  expect_equal(
    .checkColorScheme(c("viridis", "default"), c("CloneSeq", "CloneCount")),
    c("viridis", "default"), ignore_attr = TRUE
  )
  expect_equal(.checkColorScheme(NULL, NULL), NULL, ignore_attr = TRUE)
})

test_that(".checkColorTitle works correctly", {
  expect_equal(.checkColorTitle("foo", "CloneSeq"), "foo"
  )
  expect_warning(
    .checkColorTitle(c("a", "b", "c"), c("CloneSeq", "CloneCount"))
  )
  expect_warning(
    .checkColorTitle(1:2, c("CloneSeq", "CloneCount")),
  )
  expect_equal(
    .checkColorTitle("title", c("CloneSeq", "CloneCount")),
    "title"
  )
  expect_equal(.checkColorTitle(NULL, "CloneSeq"), NULL)
  expect_equal(.checkColorTitle(NULL, c("CloneSeq", "CloneCount")), NULL)
  expect_equal(.checkColorTitle(NULL, NULL), NULL)
})

test_that(".checkNodeSizeLimits works correctly", {
  expect_warning(
    .checkNodeSizeLimits(1)
  )
  expect_warning(
    .checkNodeSizeLimits(c(-1, 1))
  )
  expect_warning(
    .checkNodeSizeLimits(c(1, -1))
  )
  expect_warning(
    .checkNodeSizeLimits(c(1, 0.5))
  )
  expect_equal(.checkNodeSizeLimits(c(2, 3)), c(2, 3), ignore_attr = TRUE)
  expect_equal(.checkNodeSizeLimits(NULL), NULL, ignore_attr = TRUE)
})

test_that(".checkPlotsAgainstLayout works", {
  expect_null(.checkPlotsAgainstLayout(net$plots))
  expect_null(.checkPlotsAgainstLayout(net0$plots))
  expect_warning(.checkPlotsAgainstLayout(net2$plots),
                 "does not match one or more of the plots contained in"
  )
})
