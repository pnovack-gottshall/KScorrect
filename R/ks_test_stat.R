#' Internal KScorrect Function.
#'
#' Internal function not intended to be called directly by users.
#'
#' Simplified and faster \code{\link{ks.test}} function that calculates just the
#' two-sided test statistic D.
#'
#' @param x a numeric vector of data values.
#' @param y a character string naming a cumulative distribution function or an
#'   actual cumulative distribution function such as pnorm. Only continuous CDFs
#'   are valid. See /code{LcKS} for accepted functions.
#' @param ... parameters of the distribution specified (as a character string)
#'   by y.
#'
#' @note Calculating the Kolmogorov-Smirnov test statistic D by itself is faster
#'   than calculating the other ouput that that function produces.
#'
#' @seealso \code{\link[stats]{ks.test}}
#'
#' @export
ks_test_stat <- function(x, y, ...) {
  x <- x[!is.na(x)]
  n <- length(x)
  y <- get(y, mode = "function", envir = parent.frame())
  diffs <- y(sort(x), ...) - (0:(n - 1)) / n
  D.obs <- max(c(diffs, 1 / n - diffs))
  return(D.obs)
}
