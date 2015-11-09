#' The Normal Mixture Distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' a univariate (one-dimensional) distribution composed of a mixture of normal
#' distributions with means equal to \code{mean}, standard deviations equal to
#' \code{sd}, and mixing proportion of the components equal to \code{pro}.
#'
#' @param x Vector of quantiles.
#' @param q Vector of quantiles.
#' @param npts Number of probability points in resampling quantile
#'   approximation. \code{Default = 100} points. See \code{details} below.
#' @param nreps Number of replicates to use in resampling approximation of
#'   quantile function. \code{Default = 5000} replicates. See \code{details}
#'   below.
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
#'   identified from the length of the \code{mean} vector. If a single \code{sd}
#'   value is provided, an equal-variance mixture model (\code{modelNames="E"}
#'   in \code{mclust}) is implemented; if multiple values are provided, a
#'   variable-variance model (\code{modelNames="V"} in \code{mclust}) is
#'   implemented. If mixing proportion \code{pro} is missing, each component is
#'   assigned equal mixing proportions, with a warning called. If supplied
#'   proportions do not sum to unity, they are rescaled to do so. If the length
#'   of supplied means, standard deviations, and mixing proportions conflict, an
#'   error is called.
#'
#'   Because mixture models are not shape and scale invariant, an analytical
#'   solution is not available to calculate a quantile function for a given
#'   combination of mixture parameters. Monte Carlo resampling methods are used
#'   in \code{qmixnorm} to approximate the quantile function. The function is
#'   approximated by drawing \code{(default =) 100} random numbers from the
#'   specified mixture distribution, repeating this \code{(default =) 5000}
#'   times, and then taking the mean of the rank-ordered values to yield
#'   expected quantiles. Although 2000 \code{nreps} replicates are generally
#'   sufficient for most uses, \code{default = 5000} provides more precise
#'   approximations at slightly increased time. See \code{examples} for
#'   confirmation that approximations are accurate, comparing the approximate
#'   quantiles from a single 'mixture' distribution to those calculated
#'   analytically for the same distribution using \code{qnorm}.
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
#' # Density and distribution functions
#' plot(seq(0, 10, .1), dmixnorm(seq(0, 10, .1), mean=mean, pro=pro, sd=sd),
#'      type="l", main="Normal mixture density")
#' plot(seq(0, 10, .1), pmixnorm(seq(0, 10, .1), mean=mean, pro=pro, sd=sd),
#'      type="l", main="Normal mixture cumulative")
#'
#' # Resampling solution for quantile function
#' npts <- 100
#' nreps <- 2000
#' exp <- qmixnorm(npts, nreps, mean=mean, pro=pro, sd=sd)
#' cbind(ppoints(npts), exp)     # Quantiles for various probabilities
#' plot(exp, type="l", main="Approximate quantile for normal mixture")
#'
#' \dontrun{
#' # Requires functions from the 'mclust' package
#' # Confirmation that approximates correct solution
#' #   (single component 'mixture' should mimic qnorm):
#' x <- rmixnorm(n=5000, mean=0, pro=1, sd=1)
#' mpar <- mclust::Mclust(x)$param
#' approx <- qmixnorm(npts, nreps, mean=mpar$mean, pro=mpar$pro,
#'      sd=sqrt(mpar$variance$sigmasq))
#' known <- qnorm(ppoints(npts), mean=mpar$mean, sd=sqrt(mpar$variance$sigmasq))
#' cor(approx, known)  # Approximately the same
#' plot(approx, main="Quantiles for (unimodal) normal")
#' lines(known)
#' legend("topleft", legend=c("known", "approximation"), pch=c(NA,1),
#'      lty=c(1, NA), bty="n")
#' }
#' @export
#' @import mclust
#' @importFrom stats ppoints
dmixnorm <- function(x, mean=NULL, sd=NULL, pro=NULL) {
  if(mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector.")
  if(any(missing(mean), missing(sd)))
    stop("'mean' and 'sd' not provided, without default.")
  mean <- as.vector(mean, mode="numeric")
  G <- length(mean)
  sd <- as.vector(sd, mode="numeric")
  lsd <- length(sd)
  if (missing(pro)) {
    pro <- rep(1/G, G)
    warning("mixing proportion 'pro' not provided. Assigned equal proportions by default.")
  }
  if(any(pro < 0L, sd < 0L))
    stop("'pro' and 'sd' must not be negative.")
  lpro <- length(pro)
  if(G < lsd | G < lpro | (lsd > 1 & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode="numeric")
  pro <- pro/sum(pro)
  if(length(sd)==1) ( modelName <- "E" ) else ( modelName <- "V" )
  parameters <- list(mean=mean, pro=pro, variance=list(sigmasq=sd^2, d=1, G=G))
  dens <- mclust::dens(data=x, modelName=modelName, parameters=parameters)
  return(dens)
}
