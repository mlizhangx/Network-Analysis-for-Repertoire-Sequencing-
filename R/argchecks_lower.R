# NAIR: Network Analysis of Immune Repertoire
# Copyright (C) 2023 Li Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.



# Conventions -------------------------------------------------------------

# Checks beginning with `.is` return TRUE or FALSE
# Checks beginning with `.MUST` raise an error or return NULL (hard checks)
# Checks beginning with `.check` return the main argument or a default value
#                                                               (soft checks)


# Generic Checks ----------------------------------------------------------
# Can be called to apply other checks
# Or used as building blocks in other checks


.stopifnot <- function(condition, name, ..., sep = " ") {
  if (!isTRUE(condition)) {
    stop(paste(sQuote(name), ..., sep = sep))
  }
}

.warnifnot <- function(condition, name, ..., sep = " ") {
  if (!isTRUE(condition)) {
    warning(paste(sQuote(name), ..., sep = sep))
  }
}

.orNull <- function(check, x, ...) {
  # wrap other checks to allow NULL argument values;
  # works with hard checks (.MUST.*) and T/F returns (.is*)
  # does not work with soft checks (.check*)
  if (!is.null(x)) { return(check(x, ...)) }
  TRUE
}

.check <- function(x, check, default, ornull = FALSE, ...,
                   nse = TRUE, dquote = FALSE
) {
  # argument `check` is a function that returns TRUE or FALSE
  # Return original value if check passed, else return default value
  if ((!ornull && is.null(x)) || (!is.null(x) && !check(x, ...))) {
    if (nse) {
      default_name <- deparse(substitute(default))
    } else {
      default_name <- default
    }
    if (dquote) { default_name <- dQuote(default_name) }
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Defaulting to ", default_name
    )
    return(default)
  }
  x
}



# Universal Properties ----------------------------------------------------


.hasLength <- function(x, len) { length(x) %in% len }
.MUST.hasLength <- function(x, len, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  emsg <- ifelse(length(len) == 1,
                 yes = "must have length",
                 no = "must have one of the following lengths:"
  )
  .stopifnot(.hasLength(x, len), name, emsg, paste(len, collapse = ", "))
}

.hasPosLength <- function(x) { length(x) > 0 }
.MUST.hasPosLength <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.hasPosLength(x), name, "must have positive length")
}

.hasLength1 <- function(x) { length(x) == 1 }
.MUST.hasLength1 <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.hasLength1(x), name, "must have length 1")
}


.hasLength2 <- function(x) { length(x) == 2 }
.MUST.hasLength2 <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.hasLength2(x), name, "must have length 2")
}


.hasElem <- .hasElement <- function(x, elem) { isTRUE(elem %in% names(x)) }
.MUST.hasElem <- .MUST.hasElement <- function(x, elem, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.hasElement(x, elem),
             name, "must contain an element named", dQuote(elem)
  )
}

.hasNA <- .hasNAs <- function(x) { sum(is.na(x)) > 0 }
.noNAs <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(!.hasNA(x), name, "must not contain NA or NaN values")
}

.isFinite <- function(x) {
  is.numeric(x) && .hasPosLength(x) && all(is.finite(x))
}
.MUST.isFinite <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isFinite(x), name, "must contain finite values")
}

# Type Checks -------------------------------------------------------------



.MUST.isLogical <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(is.logical(x), name, "must be of type", dQuote("logical"))
  .noNAs(x, name)
}

.MUST.isChar <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(is.character(x), name, "must be of type", dQuote("character"))
}

.MUST.isNumeric <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(is.numeric(x), name, "must be of type", dQuote("numeric"))
}

.isCharOrNumeric <- function(x) { is.character(x) || is.numeric(x) }
.MUST.isCharOrNumeric <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isCharOrNumeric(x),
    name, "must be of type", dQuote("character"), "or", dQuote("numeric")
  )
}

.isCharOrLogical <- function(x) { is.character(x) || is.logical(x) }
.MUST.isCharOrLogical <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isCharOrLogical(x),
    name, "must be of type", dQuote("character"), "or", dQuote("logical")
  )
  if (is.logical(x)) { .noNAs(x, name) }
}



# Scalar Types ------------------------------------------------------------

.isTF <- function(x) { is.logical(x) && .hasLength1(x) && !is.na(x) }
.MUST.isTF <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isTF(x),
             name, "must evaluate to", dQuote("TRUE"), "or", dQuote("FALSE")
  )
}
.checkTF <- function(x, default, or_auto = FALSE) {
  if (or_auto && isTRUE(x == "auto")) {
    return(x)
  }
  if (!.isTF(x)) {
    warning(
      "value for ", sQuote(deparse(substitute(x))), " is invalid. ",
      "Defaulting to ", deparse(substitute(default))
    )
    return(default)
  }
  x
}

.isTFOrAuto <- function(x) {
  .isTF(x) || (is.character(x) && .hasLength1(x) && x == "auto")
}
.MUST.isTFOrAuto <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isTFOrAuto(x),
    name, "must be", dQuote("TRUE"), ",", dQuote("FALSE"), "or", dQuote("auto")
  )
}

.isNumericScalar <- function(x) { .hasLength1(x) && .isFinite(x) }
.MUST.isNumericScalar <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isNumericScalar(x), name, "must be a finite numeric scalar")
}

.isNonneg <- function(x) { .isNumericScalar(x) && isTRUE(x >= 0) }
.MUST.isNonneg <- function(x, name = NULL) {
  .stopifnot(.isNonneg(x), name, "must be a finite and nonnegative scalar")
}

.isPos <- function(x) { .isNumericScalar(x) && isTRUE(x > 0) }
.MUST.isPos <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isPos(x), name, "must be a finite and strictly positive scalar")
}

.isInt <- function(x) { .isNumericScalar(x) && isTRUE(x %% 1 == 0) }
.MUST.isInt <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isInt(x), name, "must be a finite integer")
}

.isNonnegInt <- function(x) { .isNonneg(x) && isTRUE(x %% 1 == 0) }
.MUST.isNonnegInt <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isNonnegInt(x),
             name, "must be a finite, nonnegative integer"
  )
}

.isPosInt <- function(x) { .isPos(x) && isTRUE(x %% 1 == 0) }
.MUST.isPosInt <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isPosInt(x),
             name, "must be a finite, strictly positive integer")
}

.isString <- function(x) { is.character(x) && .hasLength1(x) && !is.na(x) }
.MUST.isString <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isString(x), name, "must be a character string")
}

.isNonemptyString <- function(x) { .isString(x) && x != "" }
.MUST.isNonemptyString <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isNonemptyString(x), name, "must be a nonempty character string")
}

.isStringOrConnection <- function(x) {
  .isString(x) || inherits(x, "connection")
}
.MUST.isStringOrConnection <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isStringOrConnection(x),
             name, "must be a character string or a connection"
  )
}

.isCharOrNumericScalar <- function(x) { .isString(x) || .isNumericScalar(x) }
.MUST.isCharOrNumericScalar <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isCharOrNumericScalar(x),
             name, "must be a character string or finite numeric scalar"
  )
}

.isStringOrInt <- function(x) { .isString(x) || .isInt(x) }
.MUST.isStringOrInt <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isStringOrInt(x),
             name, "must be a character string or finite integer"
  )
}

.isStringOrPosInt <- function(x) { .isString(x) || .isPosInt(x) }
.MUST.isStringOrPosInt <- function(x, name = NULL) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isStringOrPosInt(x),
    name, "must be a character string or a finite, strictly positive integer"
  )
}


# Vector Types ------------------------------------------------------------


.isLogicalVector <- function(x) {
  is.vector(x, mode = "logical") && .hasPosLength(x)
}

.isNumericVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x, mode = "numeric") }
  is.vector(x, mode = "numeric") && .hasPosLength(x)
}

.isIntegerVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x, mode = "numeric") }
  .isNumericVector(x) && isTRUE(sum(x %% 1) == 0)
}
.MUST.isIntegerVector <- function(x, name = NULL, factor_ok = FALSE) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isIntegerVector(x, factor_ok),
             name, "must be a nonempty integer-valued vector"
  )
}

.isNonnegIntegerVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x, mode = "numeric") }
  .isNumericVector(x) && isTRUE(sum(x %% 1) == 0) && all(x >= 0)
}

.isPosIntegerVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x, mode = "numeric") }
  .isNumericVector(x) && isTRUE(sum(x %% 1) == 0) && all(x > 0)
}
.MUST.isPosIntegerVector <- function(x, name = NULL, factor_ok = FALSE) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isPosIntegerVector(x, factor_ok),
             name, "must be a nonempty vector containing",
             "strictly positive integer values"
  )
}

.isCharVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x, mode = "character") }
  is.vector(x, mode = "character") && .hasPosLength(x)
}
.MUST.isCharVector <- function(x, name = NULL, factor_ok = FALSE) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(.isCharVector(x, factor_ok), name,
             "must be a nonempty character vector"
  )
}

.isCharOrNumericVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x) }
  .hasPosLength(x) && any(is.vector(x, mode = "character"),
                          is.vector(x, mode = "numeric")
  )
}
.MUST.isCharOrNumericVector <- function(x, name = NULL, factor_ok = FALSE) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isCharOrNumericVector(x, factor_ok),
    name, "must be a nonempty vector of type",
    dQuote("character"), "or", dQuote("numeric")
  )
}

.isCharOrIntegerVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x) }
  .isCharVector(x) || .isIntegerVector(x)
}
.MUST.isCharOrIntegerVector <- function(x, name = NULL, factor_ok = FALSE) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isCharOrIntegerVector(x, factor_ok),
    name, "must be a nonempty character vector or integer-valued vector"
  )
}


.isCharOrPosIntegerVector <- function(x, factor_ok = FALSE) {
  if (factor_ok && is.factor(x)) { x <- as.vector(x) }
  .isCharVector(x) || .isPosIntegerVector(x)
}
.MUST.isCharOrPosIntegerVector <- function(x, name = NULL, factor_ok = FALSE) {
  if (is.null(name)) { name <- deparse(substitute(x)) }
  .stopifnot(
    .isCharOrPosIntegerVector(x, factor_ok),
    name, "must be a nonempty character vector or a numeric vector",
    "with strictly positive integer values"
  )
}