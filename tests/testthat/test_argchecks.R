library(NAIR)

test_that(".nonNull works correctly", {
  expect_error(
    .nonNull(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_identical(.nonNull(3, "argument"), NULL)
  expect_identical(.nonNull("foo", "argument"), NULL)
  expect_identical(.nonNull(NA, "argument"), NULL)
  expect_identical(.nonNull(logical(0), "argument"), NULL)
  expect_identical(.nonNull(numeric(0), "argument"), NULL)
  expect_identical(.nonNull(character(0), "argument"), NULL)
  expect_identical(.nonNull(NaN, "argument"), NULL)
  expect_identical(.nonNull(Inf, "argument"), NULL)
  expect_identical(.nonNull(TRUE, "argument"), NULL)
})

test_that(".noNAs works correctly", {
  expect_error(
    .noNAs(NA, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .noNAs(NaN, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .noNAs(c(3, NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .noNAs(c(3, NaN), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_identical(.noNAs(3, "argument"), NULL)
  expect_identical(.noNAs("foo", "argument"), NULL)
  expect_identical(.noNAs(NULL, "argument"), NULL)
  expect_identical(.noNAs(logical(0), "argument"), NULL)
  expect_identical(.noNAs(numeric(0), "argument"), NULL)
  expect_identical(.noNAs(character(0), "argument"), NULL)
  expect_identical(.noNAs(Inf, "argument"), NULL)
  expect_identical(.noNAs(TRUE, "argument"), NULL)
})

test_that(".hasLength1 works correctly", {
  expect_error(
    .hasLength1(NULL, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength1(logical(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength1(numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength1(character(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength1(c("fee", "fie"), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength1(diag(2), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength1(1:3, "argument"),
    "argument must have length 1"
  )
  expect_identical(.hasLength1(3, "argument"), NULL)
  expect_identical(.hasLength1("foo", "argument"), NULL)
  expect_identical(.hasLength1(NA, "argument"), NULL)
  expect_identical(.hasLength1(NaN, "argument"), NULL)
  expect_identical(.hasLength1(Inf, "argument"), NULL)
  expect_identical(.hasLength1(TRUE, "argument"), NULL)
})

test_that(".hasLength2 works correctly", {
  expect_error(
    .hasLength2(NULL, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(logical(0), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(numeric(0), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(character(0), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(diag(2), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(1:3, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(3, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2("foo", "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(NA, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(NaN, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(Inf, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength2(TRUE, "argument"),
    "argument must have length 2"
  )
  expect_identical(.hasLength2(c("fee", "fie"), "argument"), NULL)
  expect_identical(.hasLength2(c(23, NaN), "argument"), NULL)
  expect_identical(.hasLength2(c(TRUE, FALSE), "argument"), NULL)
  expect_identical(.hasLength2(numeric(2), "argument"), NULL)
})

test_that(".hasLength works correctly", {
  expect_error(
    .hasLength(2, NULL, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, logical(0), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, numeric(0), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, character(0), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, diag(2), "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, 1:3, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, 3, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, "foo", "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, NA, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, NaN, "argument"),
    "argument must have length 2"
  )
  expect_error(
    .hasLength(2, Inf, "argument"),
    "argument must have length 2"
  )
  expect_identical(.hasLength(2, c("fee", "fie"), "argument"), NULL)
  expect_identical(.hasLength(2, c(23, NaN), "argument"), NULL)
  expect_identical(.hasLength(2, c(TRUE, FALSE), "argument"), NULL)
  expect_identical(.hasLength(2, numeric(2), "argument"), NULL)
  expect_error(
    .hasLength(1, NULL, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength(1, logical(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength(1, numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength(1, character(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength(1, c("fee", "fie"), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength(1, diag(2), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .hasLength(1, 1:3, "argument"),
    "argument must have length 1"
  )
  expect_identical(.hasLength(1, 3, "argument"), NULL)
  expect_identical(.hasLength(1, "foo", "argument"), NULL)
  expect_identical(.hasLength(1, NA, "argument"), NULL)
  expect_identical(.hasLength(1, NaN, "argument"), NULL)
  expect_identical(.hasLength(1, Inf, "argument"), NULL)
  expect_error(
    .hasLength(3, 1:4, "argument"),
    "argument must have length 3"
  )
  expect_identical(.hasLength(3, 1:3, "argument"), NULL)
})

test_that(".orNull works correctly", {
  expect_error(
    .orNull(.noNAs, NA, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .orNull(.hasLength1, 1:3, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .orNull(.isNumeric, "foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .orNull(.isNumeric, "foo", "argument"),
    "argument must be of type numeric"
  )
  expect_identical(.orNull(.noNAs, NULL, "argument"), NULL)
  expect_identical(.orNull(.hasLength1, NULL, "argument"), NULL)
  expect_identical(.orNull(.isNumeric, NULL, "argument"), NULL)
})

test_that(".isLogical works correctly", {
  expect_error(
    .isLogical(NaN, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isLogical(3, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isLogical(1, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isLogical(0, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isLogical("TRUE", "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isLogical(NULL, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isLogical(NA, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_identical(.isLogical(TRUE, "argument"), NULL)
  expect_identical(.isLogical(FALSE, "argument"), NULL)
  expect_identical(.isLogical(logical(0), "argument"), NULL)
  expect_identical(.isLogical(c(TRUE, FALSE), "argument"), NULL)
})

test_that(".isChar works correctly", {
  expect_error(
    .isChar(TRUE, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isChar(3, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isChar(1, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isChar(0, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isChar(NaN, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isChar(NA, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isChar(NULL, "argument"),
    "argument must be of type character"
  )
  expect_identical(.isChar("foo", "argument"), NULL)
  expect_identical(.isChar(character(0), "argument"), NULL)
  expect_identical(.isChar(c("foo", "bar"), "argument"), NULL)
})

test_that(".isNumeric works correctly", {

  expect_error(
    .isNumeric(NULL, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNumeric(NA, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNumeric("foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNumeric(TRUE, "argument"),
    "argument must be of type numeric"
  )
  expect_identical(.isNumeric(3, "argument"), NULL)
  expect_identical(.isNumeric(numeric(0), "argument"), NULL)
  expect_identical(.isNumeric(NaN, "argument"), NULL)
  expect_identical(.isNumeric(Inf, "argument"), NULL)
  expect_identical(.isNumeric(1:4, "argument"), NULL)
})

test_that(".isCharOrNumeric works correctly", {
  expect_error(
    .isCharOrNumeric(TRUE, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isCharOrNumeric(NULL, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isCharOrNumeric(NA, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isCharOrNumeric(list(c("fee", "fie"), 1:2), "argument"),
    "argument must be of type character or numeric"
  )
  expect_identical(.isCharOrNumeric("foo", "argument"), NULL)
  expect_identical(.isCharOrNumeric(character(0), "argument"), NULL)
  expect_identical(.isCharOrNumeric(3, "argument"), NULL)
  expect_identical(.isCharOrNumeric(numeric(0), "argument"), NULL)
  expect_identical(.isCharOrNumeric(NaN, "argument"), NULL)
  expect_identical(.isCharOrNumeric(Inf, "argument"), NULL)
  expect_identical(.isCharOrNumeric(1:4, "argument"), NULL)
  expect_identical(.isCharOrNumeric(c("foo", "bar"), "argument"), NULL)
})

test_that(".isCharOrLogical works correctly", {
  expect_error(
    .isCharOrLogical(NULL, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(list(c("fee", "fie"), 1:2), "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(3, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(1, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(0, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(NaN, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(Inf, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isCharOrLogical(NA, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_identical(.isCharOrLogical("foo", "argument"), NULL)
  expect_identical(.isCharOrLogical(character(0), "argument"), NULL)
  expect_identical(.isCharOrLogical(c("foo", "bar"), "argument"), NULL)
  expect_identical(.isCharOrLogical(TRUE, "argument"), NULL)
  expect_identical(.isCharOrLogical(FALSE, "argument"), NULL)
  expect_identical(.isCharOrLogical(logical(0), "argument"), NULL)
  expect_identical(.isCharOrLogical(c(TRUE, FALSE), "argument"), NULL)
})

test_that(".hasPosLength works correctly", {
  expect_error(
    .hasPosLength(NULL, "argument"),
    "argument must have positive length"
  )
  expect_error(
    .hasPosLength(logical(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .hasPosLength(numeric(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .hasPosLength(character(0), "argument"),
    "argument must have positive length"
  )
  expect_identical(.hasPosLength(3, "argument"), NULL)
  expect_identical(.hasPosLength("foo", "argument"), NULL)
  expect_identical(.hasPosLength(NA, "argument"), NULL)
  expect_identical(.hasPosLength(NaN, "argument"), NULL)
  expect_identical(.hasPosLength(Inf, "argument"), NULL)
  expect_identical(.hasPosLength(TRUE, "argument"), NULL)
  expect_identical(.hasPosLength(c("fee", "fie"), "argument"), NULL)
  expect_identical(.hasPosLength(c(23, NaN), "argument"), NULL)
  expect_identical(.hasPosLength(c(TRUE, FALSE), "argument"), NULL)
  expect_identical(.hasPosLength(numeric(2), "argument"), NULL)
  expect_identical(.hasPosLength(1:4, "argument"), NULL)
  expect_identical(.hasPosLength(3, "argument"), NULL)
  expect_identical(.hasPosLength(diag(2), "argument"), NULL)
})

test_that(".isTF works correctly", {
  expect_error(
    .isTF(NA, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isTF(NaN, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isTF(3, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isTF(1, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isTF(0, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isTF("TRUE", "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isTF(logical(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isTF(NULL, "argument"),
    "argument must be of type logical"
  )
  expect_error(
    .isTF(c(TRUE, FALSE), "argument"),
    "argument must have length 1"
  )
  expect_identical(.isTF(TRUE, "argument"), NULL)
  expect_identical(.isTF(FALSE, "argument"), NULL)
})

test_that(".isTFOrAuto works correctly", {
  expect_error(
    .isTFOrAuto(NA, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isTFOrAuto(NaN, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isTFOrAuto(3, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isTFOrAuto(1, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isTFOrAuto(0, "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isTFOrAuto("TRUE", "argument"),
    "argument must be TRUE, FALSE, or \"auto\""
  )
  expect_error(
    .isTFOrAuto("AUTO", "argument"),
    "argument must be TRUE, FALSE, or \"auto\""
  )
  expect_error(
    .isTFOrAuto(logical(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isTFOrAuto(character(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isTFOrAuto(numeric(0), "argument"),
    "argument must be of type character or logical"
  )
  expect_error(
    .isTFOrAuto(c(TRUE, FALSE), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isTFOrAuto(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_identical(.isTFOrAuto(TRUE, "argument"), NULL)
  expect_identical(.isTFOrAuto(FALSE, "argument"), NULL)
  expect_identical(.isTFOrAuto("auto", "argument"), NULL)
})

test_that(".isFinite works correctly", {

  expect_error(
    .isFinite(NULL, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isFinite(NA, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isFinite("foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isFinite(TRUE, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isFinite(character(0), "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isFinite(numeric(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .isFinite(NaN, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isFinite(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isFinite(c(2, Inf), "argument"),
    "argument must contain finite values"
  )
  expect_identical(.isFinite(3, "argument"), NULL)
  expect_identical(.isFinite(1:4, "argument"), NULL)
})

test_that(".isNonneg works correctly", {

  expect_error(
    .isNonneg(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isNonneg(NA, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNonneg("foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNonneg(TRUE, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNonneg(character(0), "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isNonneg(numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isNonneg(1:2, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isNonneg(NaN, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isNonneg(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isNonneg(-1, "argument"),
    "argument must be nonnegative"
  )
  expect_error(
    .isNonneg(-1.5, "argument"),
    "argument must be nonnegative"
  )
  expect_identical(.isNonneg(3, "argument"), NULL)
  expect_identical(.isNonneg(0, "argument"), NULL)
})

test_that(".isPos works correctly", {

  expect_error(
    .isPos(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isPos(NA, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPos("foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPos(TRUE, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPos(character(0), "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPos(numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isPos(1:2, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isPos(NaN, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isPos(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isPos(-1, "argument"),
    "argument must be strictly positive"
  )
  expect_error(
    .isPos(0, "argument"),
    "argument must be strictly positive"
  )
  expect_error(
    .isPos(-1.5, "argument"),
    "argument must be strictly positive"
  )
  expect_identical(.isPos(3, "argument"), NULL)
})

test_that(".isInt works correctly", {

  expect_error(
    .isInt(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isInt(NA, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isInt("foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isInt(TRUE, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isInt(character(0), "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isInt(numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isInt(1:2, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isInt(NaN, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isInt(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isInt(-1.5, "argument"),
    "argument must be integer-valued"
  )
  expect_error(
    .isInt(pi, "argument"),
    "argument must be integer-valued"
  )
  expect_identical(.isInt(3, "argument"), NULL)
  expect_identical(.isInt(3.0, "argument"), NULL)
})

test_that(".isPosInt works correctly", {

  expect_error(
    .isPosInt(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isPosInt(NA, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPosInt("foo", "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPosInt(TRUE, "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPosInt(character(0), "argument"),
    "argument must be of type numeric"
  )
  expect_error(
    .isPosInt(numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isPosInt(1:2, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isPosInt(NaN, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isPosInt(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isPosInt(1.5, "argument"),
    "argument must be integer-valued"
  )
  expect_error(
    .isPosInt(-1.5, "argument"),
    "argument must be integer-valued"
  )
  expect_error(
    .isPosInt(-1, "argument"),
    "argument must be strictly positive"
  )
  expect_error(
    .isPosInt(0, "argument"),
    "argument must be strictly positive"
  )
  expect_identical(.isPosInt(3, "argument"), NULL)
  expect_identical(.isPosInt(3.0, "argument"), NULL)
})

test_that(".isString works correctly", {
  expect_error(
    .isString(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isString(TRUE, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isString(3, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isString(1, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isString(0, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isString(NaN, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isString(NA, "argument"),
    "argument must be of type character"
  )
  expect_error(
    .isString(c("foo", "bar"), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isString(character(0), "argument"),
    "argument must have length 1"
  )
  expect_identical(.isString("foo", "argument"), NULL)
})

test_that(".isCharOrNumericScalar works correctly", {
  expect_error(
    .isCharOrNumericScalar(NA, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isCharOrNumericScalar(TRUE, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isCharOrNumericScalar(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isCharOrNumericScalar(list(c("fee", "fie"), 1:2), "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isCharOrNumericScalar(1:2, "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isCharOrNumericScalar(c("fee", "fie"), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isCharOrNumericScalar(character(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isCharOrNumericScalar(numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isCharOrNumericScalar(NaN, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isCharOrNumericScalar(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_identical(.isCharOrNumericScalar("foo", "argument"), NULL)
  expect_identical(.isCharOrNumericScalar(3, "argument"), NULL)
})

test_that(".isCharVector works correctly", {
  expect_error(
    .isCharVector(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isCharVector(TRUE, "argument"),
    "argument must be a character vector"
  )
  expect_error(
    .isCharVector(3, "argument"),
    "argument must be a character vector"
  )
  expect_error(
    .isCharVector(NaN, "argument"),
    "argument must be a character vector"
  )
  expect_error(
    .isCharVector(NA, "argument"),
    "argument must be a character vector"
  )
  expect_error(
    .isCharVector(character(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .isCharVector(c("foo", NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_identical(.isCharVector("foo", "argument"), NULL)
  expect_identical(.isCharVector(c("foo", "bar"), "argument"), NULL)
})

test_that(".isCharOrNumericVector works correctly", {
  expect_error(
    .isCharOrNumericVector(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isCharOrNumericVector(TRUE, "argument"),
    "argument must be a character or numeric vector"
  )
  expect_error(
    .isCharOrNumericVector(NA, "argument"),
    "argument must be a character or numeric vector"
  )
  expect_error(
    .isCharOrNumericVector(list(c("fee", "fie"), 1:2), "argument"),
    "argument must be a character or numeric vector"
  )
  expect_error(
    .isCharOrNumericVector(character(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .isCharOrNumericVector(numeric(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .isCharOrNumericVector(c(3, NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isCharOrNumericVector(c("foo", NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isCharOrNumericVector(NaN, "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isCharOrNumericVector(c(2.1, NaN), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isCharOrNumericVector(c(2.1, NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isCharOrNumericVector(Inf, "argument"),
    "argument must contain finite values"
  )
  expect_error(
    .isCharOrNumericVector(c(2.1, Inf), "argument"),
    "argument must contain finite values"
  )
  expect_identical(.isCharOrNumericVector("foo", "argument"), NULL)
  expect_identical(.isCharOrNumericVector(c("foo", "bar"), "argument"), NULL)
  expect_identical(.isCharOrNumericVector(3, "argument"), NULL)
  expect_identical(.isCharOrNumericVector(1:4, "argument"), NULL)


})

test_that(".hasAtLeastTwoRows works correctly", {
  foo <- data.frame(a = 1, b = 1, c = 2)
  foo2 <- rbind(foo, foo)
  expect_error(
    .hasAtLeastTwoRows(foo),
    "need at least two data rows"
  )
  expect_identical(.hasAtLeastTwoRows(foo2), NULL)
})

test_that(".isDataCol and similar checks work correctly", {
  dat <- simulateToyData(sample_size = 5)
  expect_identical(.isDataCol(dat, "CloneSeq", "argument"), NULL)
  expect_identical(.isDataCol(dat, 1, "argument"), NULL)
  expect_identical(.isDataColOrNull(dat, "CloneSeq", "argument"), NULL)
  expect_identical(.isDataColOrNull(dat, NULL, "argument"), NULL)
  expect_identical(.isDataCols(dat, "CloneSeq", "argument"), NULL)
  expect_identical(
    .isDataCols(dat, c("CloneSeq", "CloneCount"), "argument"), NULL
  )
  expect_identical(
    .isDataColsOrNull(dat, c("CloneSeq", "CloneCount"), "argument"), NULL
  )
  expect_identical(.isDataColsOrNull(dat, NULL, "argument"), NULL)
  expect_error(
    .isDataCol(dat, NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isDataCol(dat, TRUE, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isDataCol(dat, NA, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isDataCol(dat, c(1, 2), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isDataCol(dat, c("CloneSeq", "CloneCount"), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isDataCol(dat, character(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isDataCol(dat, numeric(0), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .isDataCol(dat, "foo", "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCol(dat, NaN, "argument"),
    "argument is of type numeric and hence must contain finite values"
  )
  expect_error(
    .isDataCol(dat, Inf, "argument"),
    "argument is of type numeric and hence must contain finite values"
  )
  expect_error(
    .isDataCol(dat, 1.1, "argument"),
    "argument is of type numeric and hence must be integer-valued"
  )
  expect_error(
    .isDataCol(dat, 10, "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCol(dat, 0, "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCol(dat, -5, "argument"),
    "argument must specify a valid column within the data"
  )

  expect_error(
    .isDataCols(dat, NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .isDataCols(dat, TRUE, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isDataCols(dat, NA, "argument"),
    "argument must be of type character or numeric"
  )
  expect_error(
    .isDataCols(dat, character(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .isDataCols(dat, numeric(0), "argument"),
    "argument must have positive length"
  )
  expect_error(
    .isDataCols(dat, "foo", "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, NaN, "argument"),
    "argument is of type numeric and hence must contain finite values"
  )
  expect_error(
    .isDataCols(dat, Inf, "argument"),
    "argument is of type numeric and hence must contain finite values"
  )
  expect_error(
    .isDataCols(dat, 1.1, "argument"),
    "argument is of type numeric and hence must be integer-valued"
  )
  expect_error(
    .isDataCols(dat, 10, "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, 0, "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, -5, "argument"),
    "argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, c("foo", NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isDataCols(dat, c(1, NA), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isDataCols(dat, c(1, NaN), "argument"),
    "argument must not contain NA/NaNs"
  )
  expect_error(
    .isDataCols(dat, c("CloneSeq", "foo"), "argument"),
    "each element of argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, c(1, Inf), "argument"),
    "each element of argument is of type numeric and hence must contain finite values"
  )
  expect_error(
    .isDataCols(dat, c(1, 1.1), "argument"),
    "each element of argument is of type numeric and hence must be integer-valued"
  )
  expect_error(
    .isDataCols(dat, c(1, 10), "argument"),
    "each element of argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, c(1, 0), "argument"),
    "each element of argument must specify a valid column within the data"
  )
  expect_error(
    .isDataCols(dat, c(1, -10), "argument"),
    "each element of argument must specify a valid column within the data"
  )

})

test_that(".isValidSeqVector works correctly", {
    expect_error(
      .isValidSeqVector(NULL),
      "specified column or vector of receptor sequences must have positive length"
    )
    expect_error(
      .isValidSeqVector(character(0)),
      "specified column or vector of receptor sequences must have positive length"
    )
    expect_error(
      .isValidSeqVector(NA),
      "specified column or vector of receptor sequences contains only NA values after being coerced to a character vector"
    )
    expect_error(
      .isValidSeqVector(c(NA, NA)),
      "specified column or vector of receptor sequences contains only NA values after being coerced to a character vector"
    )
  expect_identical(.isValidSeqVector("foo"), NULL)
  expect_identical(.isValidSeqVector(1), NULL)
  expect_identical(.isValidSeqVector(c("foo", "bar")), NULL)
  expect_identical(.isValidSeqVector(c("foo", NA)), NULL)
  expect_identical(.isValidSeqVector(1:2), NULL)
  expect_identical(.isValidSeqVector(c(1, NA, NaN, Inf)), NULL)
  expect_identical(.isValidSeqVector(as.factor(c("foo", "bar", "bar"))), NULL)
})

test_that(".isSeqCol works correctly", {
  dat <- simulateToyData(sample_size = 5)
  expect_error(
    .isSeqCol(dat, NULL),
    "seq_col is required but value is NULL"
  )
  expect_error(
    .isSeqCol(dat, c("CloneSeq", "CloneCount", "CloneFreq")),
    "seq_col must have length 1 or 2"
  )
  expect_error(
    .isSeqCol(dat, character(0)),
    "seq_col must have length 1 or 2"
  )
  expect_error(
    .isSeqCol(dat, numeric(0)),
    "seq_col must have length 1 or 2"
  )
  expect_error(
    .isSeqCol(dat, "foo"),
    "seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, 1.1),
    "seq_col is of type numeric and hence must be integer-valued"
  )
  expect_error(
    .isSeqCol(dat, 0),
    "seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, 10),
    "seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, -5),
    "seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, c("CloneSeq", "foo")),
    "each element of seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, c(-1, -5)),
    "each element of seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, c(1, -5)),
    "each element of seq_col must specify a valid column within the data"
  )
  expect_error(
    .isSeqCol(dat, c(1, 1.1)),
    "each element of seq_col is of type numeric and hence must be integer-valued"
  )
  expect_identical(.isSeqCol(dat, "CloneSeq"), NULL)
  expect_identical(.isSeqCol(dat, "CloneCount"), NULL)
  expect_identical(.isSeqCol(dat, c("CloneSeq", "CloneCount")), NULL)
})

test_that(".hasElement works correctly", {
  foo <- c("first" = 1, "second" = 2, "third" = 3)
  expect_error(
    .hasElement(foo, "argument", "fourth"),
    "argument does not contain an element named fourth"
  )
  expect_error(
    .hasElement(NULL, "argument", "first"),
    "argument does not contain an element named first"
  )
  expect_error(
    .hasElement(1:3, "argument", "first"),
    "argument does not contain an element named first"
  )
  expect_error(
    .hasElement(c("first", "second", "third"), "argument", "first"),
    "argument does not contain an element named first"
  )
  expect_identical(.hasElement(foo, "argument", "first"), NULL)
})

test_that("network output type checks work correctly", {
  dat <- simulateToyData()
  net <- buildRepSeqNetwork(
    dat, "CloneSeq", cluster_stats = TRUE, plots = FALSE, output_dir = NULL
  )
  net$plots <- list("first" = ggraph::ggraph(net$igraph))
  expect_error(
    .isIgraph(net$node_data, "argument"),
    "argument must be of class igraph"
  )
  expect_error(
    .isGgraph(net$plots, "argument"),
    "argument must be of class ggraph"
  )
  expect_error(
    .isDataFrame(net$plots, "argument"),
    "argument must be a data frame"
  )
  expect_error(
    .isList(net$adjacency_matrix, "argument"),
    "argument must be a list"
  )
  expect_error(
    .isAdjacencyMatrix(net, "argument"),
    "argument must be a matrix or dgCMatrix"
  )
  expect_error(
    .isAdjacencyMatrix(matrix(0, nrow = 1, ncol = 2), "argument"),
    "argument must have the same row and column dimensions"
  )
  expect_error(
    .isAdjacencyMatrix(matrix(2, nrow = 2, ncol = 2), "argument"),
    "argument contains values other than 0 or 1"
  )
  expect_error(
    .isBaseNetworkOutput(net$node_data, "argument"),
    "argument does not contain an element named node_data"
  )
  expect_error(
    .checkIgraphAgainstData(1, diag(2)),
    "number of nodes in igraph does not match number of rows in data"
  )
  expect_error(
    .checkIgraphAgainstMatrix(1, diag(2)),
    "number of nodes in igraph does not match dimensions of adjacency matrix"
  )
  expect_error(
    .checkDataAgainstMatrix(diag(3), diag(2)),
    "number of data rows does not match dimensions of the adjacency matrix"
  )
  expect_identical(
    .isAdjacencyMatrix(matrix(0, nrow = 2, ncol = 2), "argument"), NULL
  )
  expect_identical(.isBaseNetworkOutput(net, "argument"), NULL)
  expect_identical(.hasNodeAndClusterData(net, "argument"), NULL)
  expect_identical(.isPlotlist(net$plots, "argument"), NULL)
})

test_that(".isDistType works correctly", {
  expect_error(
    .isDistType("foo"),
    "Invalid option for dist_type argument"
  )
  expect_identical(.isDistType("lev"), NULL)
})

test_that(".isInputType works correctly", {
  valid_input_types <- c("csv", "table", "tsv", "txt", "rds", "rda")
  for (i in 1:length(valid_input_types)) {
    expect_identical(.isInputType(valid_input_types[[i]]), NULL)
  }
  valid_input_types <- paste(valid_input_types, collapse = ", ")
  expect_error(
    .isInputType("foo"),
    paste("input_type must be one of:", valid_input_types)
  )
})

test_that(".isOutputType works correctly", {
  expect_warning(
    .isOutputType("csv"),
    "output_type is invalid. Defaulting to rda"
  )
  expect_warning(
    .isOutputType("tsv", "findPublicClusters"),
    "output_type is invalid. Defaulting to rds"
  )
  expect_warning(
    .isOutputType("individual", "findPublicClusters"),
    "output_type is invalid. Defaulting to rds"
  )
  expect_warning(
    .isOutputType("individual", "findAssociatedClones"),
    "output_type is invalid. Defaulting to csv"
  )
  expect_warning(
    .isOutputType("table", "findAssociatedClones"),
    "output_type is invalid. Defaulting to csv"
  )
  expect_warning(
    .isOutputType("individual", "generic"),
    "output_type is invalid. Defaulting to rda"
  )
  expect_identical(
    .isOutputType("individual"), NULL
  )
  expect_identical(
    .isOutputType("rds", "findPublicClusters"), NULL
  )
  expect_identical(
    .isOutputType("tsv", "findAssociatedClones"), NULL
  )
  expect_identical(
    .isOutputType("table", "generic"), NULL
  )
})

test_that(".checkColorNodesBy works correctly", {
  dat <- simulateToyData(sample_size = 5)
  expect_error(
    .checkColorNodesBy("foo", dat),
    "color_nodes_by specifies one or more variables not present in data or among the node-level network properties to be computed"
  )
  expect_error(
    .checkColorNodesBy("degree", dat),
    "color_nodes_by specifies one or more variables not present in data or among the node-level network properties to be computed"
  )
  expect_error(
    .checkColorNodesBy(c("CloneSeq", "degree"), dat),
    "color_nodes_by specifies one or more variables not present in data or among the node-level network properties to be computed"
  )
  expect_identical(
    .checkColorNodesBy(c("CloneSeq", "CloneCount"), dat),
    NULL
  )
  expect_identical(
    .checkColorNodesBy("CloneSeq", dat),
    NULL
  )
  expect_identical(
    .checkColorNodesBy("auto", dat),
    NULL
  )
  expect_identical(
    .checkColorNodesBy("degree", dat, node_stats = TRUE),
    NULL
  )
  expect_identical(
    .checkColorNodesBy(c("CloneSeq", "degree"), dat, node_stats = TRUE),
    NULL
  )
  expect_identical(
    .checkColorNodesBy("foo", dat, plots = FALSE),
    NULL
  )
  expect_identical(
    .checkColorNodesBy(NULL, dat),
    NULL
  )
})

test_that(".checkColorScheme works correctly", {
  expect_error(
    .checkColorScheme("foo", "CloneSeq"),
    "color_scheme contains one or more values which are not supported"
  )
  expect_error(
    .checkColorScheme(c("viridis", "foo"), c("CloneSeq", "CloneCount")),
    "color_scheme contains one or more values which are not supported"
  )
  expect_error(
    .checkColorScheme(c("foo", "bar", "baz"), c("CloneSeq", "CloneCount")),
    "color_scheme must have length 1 or the same length as color_nodes_by"
  )
  expect_null(
    .checkColorScheme("viridis", c("CloneSeq", "CloneCount"))
  )
  expect_null(
    .checkColorScheme(c("viridis", "default"), c("CloneSeq", "CloneCount"))
  )
  expect_null(
    .checkColorScheme(NULL, NULL)
  )
  expect_null(
    .checkColorScheme("foo", "CloneSeq", plots = FALSE)
  )
})

test_that(".checkSizeNodesBy works correctly", {
  dat <- simulateToyData(sample_size = 5)
  expect_error(
    .checkSizeNodesBy(1:2, dat),
    "size_nodes_by is non-null and hence must have length 1"
  )
  expect_error(
    .checkSizeNodesBy("foo", dat),
    "size_nodes_by is of type character and hence must specify a valid column within the data"
  )
  expect_error(
    .checkSizeNodesBy(0, dat),
    "size_nodes_by is of type numeric and hence must be strictly positive"
  )
  expect_null(
    .checkSizeNodesBy("CloneSeq", dat)
  )
  expect_null(
    .checkSizeNodesBy(1.5, dat)
  )
  expect_null(
    .checkSizeNodesBy(NULL, dat)
  )
})

test_that(".checkNodeSizeLimits works correctly", {
  expect_error(
    .checkNodeSizeLimits(1),
    "node_size_limits is non-null and hence must have length 2"
  )
  expect_error(
    .checkNodeSizeLimits(c(-1, 1)),
    "values for node_size_limits must be strictly positive"
  )
  expect_error(
    .checkNodeSizeLimits(c(1, -1)),
    "values for node_size_limits must be strictly positive"
  )
  expect_error(
    .checkNodeSizeLimits(c(1, 0.5)),
    "first entry of node_size_limits cannot be greater than the second entry"
  )
  expect_null(
    .checkNodeSizeLimits(c(2, 3))
  )
  expect_null(
    .checkNodeSizeLimits(NULL)
  )
})

test_that(".checkStatsToInclude works correctly", {
  expect_error(
    .checkStatsToInclude(NULL),
    "stats_to_include is required but value is NULL"
  )
  expect_error(
    .checkStatsToInclude("foo"),
    "value for stats_to_include does not match required format. See help file for chooseNodeStats()"
  )
  expect_error(
    .checkStatsToInclude(c("degree" = TRUE, "cluster_id" = TRUE)),
    "value for stats_to_include does not match required format. See help file for chooseNodeStats()"
  )
  expect_null(
    .checkStatsToInclude(chooseNodeStats())
  )
  expect_null(
    .checkStatsToInclude(exclusiveNodeStats())
  )
  expect_null(
    .checkStatsToInclude("all")
  )
  expect_null(
    .checkStatsToInclude("cluster_id_only")
  )
})

test_that(".checkClusterFun works correctly", {
  expect_error(
    .checkClusterFun(NULL, "argument"),
    "argument is required but value is NULL"
  )
  expect_error(
    .checkClusterFun(diag(5), "argument"),
    "argument must have length 1"
  )
  expect_error(
    .checkClusterFun(log, "argument"),
    "argument must be a valid clustering algorithm. See help topic \"clustering_algorithms\""
  )
  expect_error(
    .checkClusterFun("log", "argument"),
    "argument must be a valid clustering algorithm. See help topic \"clustering_algorithms\""
  )
  expect_error(
    .checkClusterFun(NA, "argument"),
    "argument must be a valid clustering algorithm. See help topic \"clustering_algorithms\""
  )
  expect_null(
    .checkClusterFun(cluster_fast_greedy, "argument")
  )
  expect_null(
    .checkClusterFun(cluster_walktrap, "argument")
  )
  expect_null(
    .checkClusterFun("cluster_fast_greedy", "argument")
  )
})
