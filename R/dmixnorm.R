#' The Normal Mixture Distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' a univariate (one-dimensional) distribution composed of a mixture of normal
#' distributions with means equal to \code{mean}, standard deviations equal to
#' \code{sd}, and mixing proportion of the components equal to \code{pro}.
#'
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param nr Number of random observations to use in quantile approximation.
#'   \code{Default = 200000} numbers. See \code{details} below.
#' @param n Number of observations.
#' @param mean Vector of means, one for each component.
#' @param sd Vector of standard deviations. If a single value is provided, an
#'   equal variance mixture model is implemented. Must not be negative.
#' @param pro Vector of mixing proportions, one for each component. If missing,
#'   an equal-proportion model is implemented, with a warning. If proportions do
#'   not sum to unity, they are rescaled to do so. Must not be negative.
#'
#' @details These functions use, modify, and wrap around those from the
#'   \code{mclust} package, especially \code{\link[mclust]{dens}} and
#'   \code{\link[mclust]{sim}}. \code{rmixnorm} is slightly faster than
#'   \code{sim} when used with univariate distributions.
#'
#'   The number of mixture components (argument \code{G} in \code{mclust}) is
#'   specified from the length of the \code{mean} vector. If a single \code{sd}
#'   value is provided, an equal-variance mixture model (\code{modelNames="E"}
#'   in \code{mclust}) is implemented; if multiple values are provided, a
#'   variable-variance model (\code{modelNames="V"} in \code{mclust}) is
#'   implemented. If mixing proportion \code{pro} is missing, all components are
#'   assigned equal mixing proportions, with a warning called. If supplied
#'   proportions do not sum to unity, they are rescaled to do so. If the length
#'   of supplied means, standard deviations, and mixing proportions conflicts,
#'   an error is called.
#'
#'   Analytical solutions are not available to calculate a quantile function for
#'   all combinations of mixture parameters. \code{qmixnorm} approximates the
#'   quantile function by drawing \code{(default) nr = 200000} random numbers
#'   from the specified mixture distribution, and using \code{stats::quantile}
#'   to produce expected quantiles for the specified probabilities \code{p}.
#'   Sensitivity analyses demontrate that using \code{[default] nr = 200000}
#'   random numbers will provide quantile values with +/- 0.01 precision 99% of
#'   the time. Using \code{[default] nr = 10000000} random numbers will provide
#'   quantile values with +/- 0.001 precision 97% of the time (at greater
#'   computational cost). See \code{examples} for confirmation that
#'   approximations are accurate, comparing the approximate quantiles from a
#'   single 'mixture' distribution to those calculated analytically for the same
#'   distribution using \code{qnorm}.
#'
#' @return \code{dmixnorm} gives the density, \code{pmixnorm} gives the
#'   distribution function, \code{qmixnorm} approximates the quantile function,
#'   and \code{rmixnorm} generates random numbers.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}, with
#'   assistance from Luca Scrucca for \code{pmixnorm}.
#' @author Steve Wang \email{scwang@@swarthmore.edu} for \code{qmixnorm}.
#'
#' @seealso \code{\link[stats]{Distributions}} for other standard distributions,
#'   and \code{mclust::\link[mclust]{dens}} and \code{\link[mclust]{sim}} for
#'   alternative density and random number functions for multivariate mixture
#'   distributions.
#'
#' @examples
#' # Mixture of two normal distributions
#' mean <- c(3, 6)
#' pro <- c(.25, .75)
#' sd <- c(.5, 1)
#' x <- rmixnorm(n=5000, mean=mean, pro=pro, sd=sd)
#' hist(x, n=20, main="random bimodal sample")
#'
#' \dontrun{
#' # Requires functions from the 'mclust' package
#' require(mclust)
#' # Confirm 'rmixnorm' above produced specified model
#' mod <- mclust::Mclust(x)
#' mod             # Best model (correctly) has two-components with unequal variances
#' mod$parameters	# and approximately same parameters as specified above
#' sd^2            # Note reports var (sigma-squared) instead of sd used above
#' }
#'
#' # Density, distribution, and quantile functions
#' plot(seq(0, 10, .1), dmixnorm(seq(0, 10, .1), mean=mean, sd=sd, pro=pro),
#'      type="l", main="Normal mixture density")
#' plot(seq(0, 10, .1), pmixnorm(seq(0, 10, .1), mean=mean, sd=sd, pro=pro),
#'      type="l", main="Normal mixture cumulative")
#' plot(seq(0, 1, .01), qmixnorm(seq(0, 1, .01), mean=mean, sd=sd, pro=pro),
#'      type="l", main="Normal mixture quantile")
#'
#' \dontrun{
#' # Requires functions from the 'mclust' package
#' # Confirmation that qmixnorm approximates correct solution
#' #   (single component 'mixture' should mimic qnorm):
#' x <- rmixnorm(n=5000, mean=0, pro=1, sd=1)
#' mpar <- mclust::Mclust(x)$param
#' approx <- qmixnorm(p=ppoints(100), mean=mpar$mean, pro=mpar$pro,
#'      sd=sqrt(mpar$variance$sigmasq))
#' known <- qnorm(p=ppoints(100), mean=mpar$mean, sd=sqrt(mpar$variance$sigmasq))
#' cor(approx, known)  # Approximately the same
#' plot(approx, main="Quantiles for (unimodal) normal")
#' lines(known)
#' legend("topleft", legend=c("known", "approximation"), pch=c(NA,1),
#'      lty=c(1, NA), bty="n")
#' }
#' @export
#' @import mclust
dmixnorm <- function(x, mean, sd, pro) {
  if(mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector.")
  if(any(missing(mean), missing(sd)))
    stop("'mean' and 'sd' not provided, without default.")
  mean <- as.vector(mean, mode="numeric")
  G <- length(mean)
  sd <- as.vector(sd, mode="numeric")
  if (missing(pro)) {
    pro <- rep(1/G, G)
    warning("mixing proportion 'pro' not provided. Assigned equal proportions by default.")
  }
  if(any(pro < 0L, sd < 0L))
    stop("'pro' and 'sd' must not be negative.")
  lpro <- length(pro)
  modelName = "V"
  lsd <- length(sd)
  if(lsd==1) {
    modelName <- "E"
    sd[seq(G)] <- sd[1]
    lsd <- length(sd)
    warning("'equal variance model' implemented. If want 'variable-variance model', specify remaining 'sd's.")
  }
  if(G < lsd | G < lpro | (lsd > 1 & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode="numeric")
  pro <- pro/sum(pro)
  parameters <- list(mean=mean, pro=pro, variance=list(sigmasq=sd^2, d=1, G=G))
  dens <- mclust::dens(data=x, modelName=modelName, parameters=parameters)
  return(dens)
}
