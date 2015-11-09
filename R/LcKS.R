#' Lilliefors-corrected Kolmogorov-Smirnoff Goodness-Of-Fit Test
#'
#' Implements the Lilliefors-corrected Kolmogorov-Smirnoff test for use in
#' goodness-of-fit tests, suitable when using sample statistics as estimates of
#' population parameters. It uses a Monte Carlo resampling algorithm to estimate
#' \emph{p}-values (and their standard error) from the resampling distribution.
#' Coded to wrap around \code{\link[stats]{ks.test}}, it is able to be used with
#' a variety of continuous distributions, including normal, lognormal,
#' univariate mixtures of normals, uniform, loguniform, exponential, gamma, and
#' Weibull distributions.
#'
#' @param x A numeric vector of data values (observed sample).
#' @param cdf Character string naming a cumulative distribution function. Case
#'   insensitive. Only continuous CDFs are valid. Allowed cdfs
#'   include:\itemize{\item \code{"pnorm"} for normal, \item \code{"pmixnorm"}
#'   for (univariate) normal mixture, \item \code{"plnorm"} for lognormal
#'   (log-normal, log normal), \item \code{"punif"} for uniform, \item
#'   \code{"plunif"} for loguniform (log-uniform, log uniform), \item
#'   \code{"pexp"} for exponential, \item \code{"pgamma"} for gamma, \item
#'   \code{"pweibull"} for Weibull.}
#' @param nreps Number of replicates to use in resampling algorithm.
#'   \code{Default = 4999} replicates. See \code{details} below. Should be a
#'   positive integer.
#' @param G Numeric vector of mixture components to consider, for mixture models
#'   only. \code{Default = 1:9} fits up to 9 components. Must contain positive
#'   integers less than 10. See \code{details} below.
#'
#'
#' @details The function builds a resampling distribution \code{D.resample} of
#'   length \code{nreps} by drawing random samples from the specified continuous
#'   distribution function \code{cdf} with parameters calculated from the
#'   provided sample \code{x}. Observed statistic \code{D} and resampled test
#'   statistics are calculated using \code{\link[stats]{ks.test}}.
#'
#'   The default \code{nreps=4999} provides generally consistent \emph{p}-values
#'   (with small standard errors). \code{nreps=1999} is sufficient for most
#'   cases, and computationally faster when dealing with more complicated
#'   distributions (such as univariate normal mixtures, gamma, and Weibull).
#'
#'   The \emph{p}-value is calculated as the number of Monte Carlo samples with
#'   D test statistics more extreme than that in the observed sample
#'   \code{D.obs}, divided by the \code{nreps} number of Monte Carlo samples. A
#'   value of 1 is added to both the numerator and denominator to allow the
#'   observed sample to be represented within the null distribution (Manly
#'   2004); this has the benefit of avoiding nonsensical \code{p.value=0.000}
#'   and accounts for the fact that the \emph{p}-value is an estimate, with an
#'   associated standard error \code{se}.
#'
#'   Sample statistics are calculated based on the specified continuous
#'   distribution, using maximum-likelihood estimates. When testing against the
#'   gamma and Weibull distributions, \code{MASS::\link[MASS]{fitdistr}} is used
#'   to calculate sample statistics using maximum likelihood optimization, with
#'   sensible starting values. Because this incorporates an optimization
#'   routine, the resampling algorithm can be slow if using large \code{nreps}
#'   or problematic samples. Warnings often occur during these optimizations,
#'   caused by difficulties estimating sample statistic standard errors. Because
#'   such SEs are not used in the Lilliefors-corrected resampling algorithm,
#'   warnings are suppressed during these optimizations.
#'
#'   Sample statistics for the (univariate) normal mixture distribution
#'   \code{\link{pmixnorm}} are calculated using package \code{mclust}, which
#'   uses BIC to identify the optimal mixture model for the sample, and the EM
#'   algorithm to estimate sample statistics for this model. The number of
#'   mixture components \code{G} (with default allowing up to 9 components),
#'   variance model (whether equal \code{E} or variable \code{V} variance), and
#'   component statistics (\code{means, sds}, and mixing proportions \code{pro})
#'   are estimated from the sample when calculating \code{D.obs} and passed
#'   internally when creating random Monte Carlo samples. It is possible that
#'   some of these samples may differ in their optimal \code{G} (for example a
#'   two-component input sample might yield a three-component random sample
#'   within the resampling distribution). This can be constrained by specifying
#'   that resampling BIC-optimizations only consider \code{G} mixture
#'   components.
#'
#'   But be aware that constraining \code{G} changes the null hypothesis. The
#'   default (\code{G=1:9}) null hypothesis is that a sample was drawn from
#'   \emph{any \code{G=1:9}-component mixture distribution}. Specifying a
#'   particular value, such as \code{G=2}, restricts the null hypothesis to
#'   particular mixture distributions with just \code{G} components, even if
#'   resampled samples might better be represented as different mixture models.
#'
#'   The \code{LcKS(cdf="pmixnorm")} test implements two control loops to avoid
#'   errors caused by this constraint and when working with problematical
#'   samples. The first loop occurs during model-selection for the observed
#'   sample \code{x}, and allows for estimation of parameters for the
#'   second-best model when those for the optimal model are not able to be
#'   calculated by the EM algorithm. A second loop occurs during the resampling
#'   algorithm, rejecting samples that can not be fit by the mixture model
#'   specified by the observed sample \code{x}. Such problematic cases are most
#'   common when the observed or resampled samples have a component(s) with very
#'   small variance (i.e., duplicate observations) or when a Monte Carlo sample
#'   can not be fit by the specified \code{G}.
#'
#' @return A list containing the following components:
#'
#'   \item{D.obs}{The value of the test statistic D for the observed sample.}
#'   \item{D.resample}{Resampling distribution of test statistics, with
#'   \code{length=nreps}. This can be used to calculate critical values; see
#'   examples.} \item{p.value}{\emph{p}-value of the test, calculated as
#'   \eqn{(\sum(D.resampled > D.obs) + 1) / (nreps + 1)}.} \item{se}{Standard
#'   error for the \emph{p}-value, calculated as \eqn{\sqrt(p.value * (1 -
#'   p.value) / nreps)}.}
#'
#' @note The Kolmogorov-Smirnoff (such as \code{ks.test}) is only valid as a
#'   goodness-of-fit test when the population parameters are known absolutely.
#'   This is typically not the case in practice. This invalidation occurs
#'   because estimating the parameters changes the null distribution of the test
#'   statistic; i.e., using the sample to estimate population parameters brings
#'   the Kolmogorov-Smirnoff test statistic D closer to the null distribution
#'   than it would be under the hypothesis where the population parameters are
#'   unknown (Gihman 1952). In other words, it is biased with a reduced Type-I
#'   error. Lilliefors (1967, 1969) provided a solution, using Monte Carlo
#'   methods to approximate the shape of the null distribution when the sample
#'   statistics are used as population parameters, and to use this null
#'   distribution as the basis for critical values. The function \code{LcKS}
#'   generalizes this solution for a range of continuous distributions.

#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}, based on
#'   code from Charles Geyer (University of Minnesota).
#'
#' @references Gihman, I. I. 1952.  Ob empiriceskoj funkcii raspredelenija
#'   slucaje grouppirovki dannych [On the empirical distribution function in the
#'   case of grouping of observations]. \emph{Doklady Akademii Nauk SSSR} 82:
#'   837-840.
#' @references Lilliefors, H. W. 1967. On the Kolmogorov-Smirnov test for
#'   normality with mean and variance unknown. \emph{Journal of the American
#'   Statistical Association} 62(318):399-402.
#' @references Lilliefors, H. W. 1969. On the Kolmogorov-Smirnov test for the
#'   exponential distribution with mean unknown. \emph{Journal of the American
#'   Statistical Association} 64(325):387-389.
#' @references Manly, B. F. J. 2004. \emph{Randomization, Bootstrap and Monte
#'   Carlo Methods in Biology}. Chapman & Hall, Cornwall, Great Britain.
#' @references Parsons, F. G., and P. H. Wirsching. 1982. A Kolmogorov-Smirnov
#'   goodness-of-fit test for the two-parameter Weibull distribution when the
#'   parameters are estimated from the data. \emph{Microelectronics Reliability}
#'   22(2):163-167.
#'
#' @seealso \code{\link[stats]{Distributions}} for standard cumulative
#'   distribution functions, \code{\link{plunif}} for the loguniform cumulative
#'   distribution function, and \code{\link{pmixnorm}} for the univariate normal
#'   mixture cumulative distribution function.
#'
#' @examples
#' x <- runif(200)
#' l <- LcKS(x, cdf="pnorm", nreps=999)
#' hist(l$D.resample)
#' abline(v = l$D.obs, lty = 2)
#' print(l, max=50)  # Just print first 50 resampled statistics
#' # Approximate p-value (usually) << 0.05
#'
#' # Confirmation uncorrected version has biased type-I error as a one-sample
#' #   test using sample statistics for parameters:
#' ks.test(x, "pnorm", mean(x), sd(x))   # p-value always larger, (usually) > 0.05
#'
#' # Confirm critical values for normal distribution are correct
#' nreps <- 9999
#' x <- rnorm(25)
#' l <- LcKS(x, "pnorm", nreps=nreps)
#' res.Ds <- sort(l$D.resample)
#' crit <- round(c(.8, .85, .9, .95, .99) * nreps, 0)
#' # Lilliefors' (1967) critical values, using improved values from
#' #   Parsons & Wirsching (1982) (for n=25):
#' # 0.141 0.148 0.157 0.172 0.201
#' round(res.Ds[crit], 3)			# Approximate critical values (converges as nreps increases)
#'
#' # Confirm critical values for exponential the same as reported by Lilliefors (1969)
#' nreps <- 9999
#' x <- rexp(25)
#' l <- LcKS(x, "pexp", nreps=nreps)
#' res.Ds <- sort(l$D.resample)
#' crit <- round(c(.8, .85, .9, .95, .99) * nreps, 0)
#' # Lilliefors' (1969) critical values (for n=25):
#' # 0.170 0.180 0.191 0.210 0.247
#' round(res.Ds[crit], 3)			# Approximate critical values (converges as nreps increases)
#'
#' \dontrun{
#' # Gamma and Weibull tests require functions from the 'MASS' package
#' # Takes time for maximum likelihood optimization of statistics
#' require(MASS)
#' x <- runif(100, min=1, max=100)
#' l <- LcKS(x, cdf="pgamma", nreps=499)
#' l$p.value
#'
#' # Confirm critical values for Weibull the same as reported by Parsons & Wirsching (1982)
#' nreps <- 9999
#' x <- rweibull(25, shape=1, scale=1)
#' l <- LcKS(x, "pweibull", nreps=nreps)
#' res.Ds <- sort(l$D.resample)
#' crit <- round(c(.8, .85, .9, .95, .99) * nreps, 0)
#' # Parsons & Wirsching (1982) critical values (for n=25):
#' # 0.141 0.148 0.157 0.172 0.201
#' round(res.Ds[crit], 3)			# Approximate critical values (converges as nreps increases)
#'
#' # Mixture test requires functions from the 'mclust' package
#' # Takes time to identify model parameters
#' require(mclust)
#' x <- rmixnorm(200, mean=c(10, 20), sd=2, pro=c(1,3))
#' l <- LcKS(x, cdf="pmixnorm", nreps=499, G=1:9)   # Default G (1:9) takes long time
#' l$p.value
#' G <- Mclust(x)$parameters$variance$G             # Optimal model has only two components
#' l <- LcKS(x, cdf="pmixnorm", nreps=499, G=G)     # Restricting to likely G saves time
#' # But note changes null hypothesis: now testing against just two-component mixture
#' l$p.value
#' }
#'
#' @export
#' @import mclust
#' @import MASS
#' @importFrom stats sd
#' @importFrom stats cor
#' @importFrom stats rnorm
#' @importFrom stats rlnorm
#' @importFrom stats runif
#' @importFrom stats rweibull
#' @importFrom stats rexp
#' @importFrom stats rgamma
#' @importFrom stats ks.test
LcKS <- function(x, cdf, nreps=4999, G=1:9) {
  if (missing(x) || length(x) == 0L || mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x)))
    stop("'x' contains missing or infinite values")
  if (nreps < 1)
    stop("'nreps' must be a positive integer.")
  if (missing(cdf) || !is.character(cdf))
    stop("'cdf' must be a character string.")
  cdf <- tolower(cdf)
  allowed <- c("pnorm", "plnorm", "punif", "plunif", "pmixnorm", "pexp", "pgamma", "pweibull")
  if(!(cdf %in% allowed))
    stop("'cdf' is not a supported distribution function.")
  if (length(G) == 1L && G == 1 && cdf == "pmixnorm")
    stop("'G' supplied not consistent with mixture model. Use 'pnorm' or 'plnorm' instead.")
  if (!identical(G, 1:9) && cdf != "pmixnorm")
    warning("'G' is ignored except when cdf='pmixnorm'.")
  if (any(G < 1 || G > 9))
    stop("'G' must be integer or vector of integers spanning 1:9.")
  n <- length(x)
  D.resample <- rep(NA, nreps)
  if(cdf=="pnorm") {
    mean.x <- mean(x)
    sd.x <- sd(x)
    D.obs <- as.vector(ks.test(x, "pnorm", mean=mean.x, sd=sd.x)$statistic)
    for (i in 1:nreps) {
      x.resample <- rnorm(n, mean=mean.x, sd=sd.x)
      D.resample[i] <- as.vector(ks.test(x.resample, "pnorm", mean=mean(x.resample), sd=sd(x.resample))$statistic)
    }
  }
  if(cdf=="plnorm") {
    if (any(x <= 0))
      stop("'x' values must be > 0 for lognormal distributions.")
    meanlog.x <- mean(log(x))
    sdlog.x <- sd(log(x))
    D.obs <- as.vector(ks.test(x, "plnorm", meanlog=meanlog.x, sdlog=sdlog.x)$statistic)
    for (i in 1:nreps) {
      x.resample <- rlnorm(n, meanlog=meanlog.x, sdlog=sdlog.x)
      D.resample[i] <- as.vector(ks.test(x.resample, "plnorm", meanlog=mean(log(x.resample)), sdlog=sd(log(x.resample)))$statistic)
    }
  }
  if(cdf=="punif") {
    min.x <- min(x)
    max.x <- max(x)
    D.obs <- as.vector(ks.test(x, "punif", min=min.x, max=max.x)$statistic)
    for (i in 1:nreps) {
      x.resample <- runif(n, min=min.x, max=max.x)
      D.resample[i] <- as.vector(ks.test(x.resample, "punif", min=min(x.resample), max=max(x.resample))$statistic)
    }
  }
  if(cdf=="plunif") {
    if (any(x <= 0))
      stop("'x' values must be > 0 for loguniform distributions.")
    min.x <- min(x)
    max.x <- max(x)
    D.obs <- as.vector(ks.test(x, "plunif", min=min.x, max=max.x)$statistic)
    for (i in 1:nreps) {
      x.resample <- rlunif(n, min=min.x, max=max.x)
      D.resample[i] <- as.vector(ks.test(x.resample, "plunif", min=min(x.resample), max=max(x.resample))$statistic)
    }
  }
  if(cdf=="pmixnorm") {
    G <- as.vector(G, mode="numeric")
    m <- mclust::mclustBIC(x, G=G)
    m <- mclust::pickBIC(m, k = sum(!is.na(m)))
    listofmod <- strsplit(names(m), ",")
    for(lom in 1:length(listofmod)) {
      mixnorm <- mclust::Mclust(x, G=listofmod[[lom]][2], modelNames=listofmod[[lom]][1], control=emControl(eps=1e-320))
      if(!all(is.na(mixnorm$parameters$pro))) break()
    }
    parameters <- mclust::Mclust(x, G=G)$parameters
    modelName <- parameters$variance$modelName
    mean.x <- parameters$mean
    sd.x <- sqrt(parameters$variance$sigmasq)
    pro.x <- parameters$pro
    D.obs <- as.vector(ks.test(x, "pmixnorm", mean=mean.x, sd=sd.x, pro=pro.x)$statistic)
    if (parameters$variance$G == 1)
      warning("Optimal mixture model for supplied sample has a single component: it is not a mixture model. Use 'pnorm' or 'plnorm' instead.")
    i <- 0
    while(i < nreps) {
      x.resample <- rmixnorm(n, mean=mean.x, pro=pro.x, sd=sd.x)
      attempt <- try(mclust::Mclust(x.resample, G=G)$parameters, silent=TRUE)
      if(inherits(attempt, "try-error")) { next }
      param.res <- attempt
      i <- i + 1
      mean.res <- param.res$mean
      sd.res <- sqrt(param.res$variance$sigmasq)
      pro.res <- param.res$pro
      D.resample[i] <- as.vector(ks.test(x.resample, "pmixnorm", mean=mean.res, sd=sd.res, pro=pro.res)$statistic)
    }
  }
  if(cdf=="pexp") {
    if (any(x <= 0))
      stop("'x' values must be > 0 for exponential distributions.")
    rate.x <- 1 / mean(x)
    D.obs <- as.vector(ks.test(x, "pexp", rate=rate.x)$statistic)
    for (i in 1:nreps) {
      x.resample <- rexp(n, rate=rate.x)
      D.resample[i] <- as.vector(ks.test(x.resample, "pexp", rate=(1/mean(x.resample)))$statistic)
    }
  }
  if(cdf=="pgamma") {
    if (any(x <= 0))
      stop("'x' values must be > 0 for gamma distributions.")
    shape.x <- mean(x)^2 / var(x)
    scale.x <- var(x) / mean(x)
    # Re-scale in case the scale is much larger than 1:
    param <- as.vector(suppressWarnings(MASS::fitdistr(x/scale.x, densfun="gamma", start=list(shape=shape.x, scale=scale.x/scale.x), control=list(maxit=25000))$estimate))
    param[2] <- param[2] * scale.x	# Re-scale back
    D.obs <- as.vector(ks.test(x, "pgamma", shape=param[1], scale=param[2])$statistic)
    for (i in 1:nreps) {
      x.resample <- rgamma(n, shape=param[1], scale=param[2])
      shape.res <- mean(x.resample)^2 / var(x.resample)
      scale.res <- var(x.resample) / mean(x.resample)
      param.res <- as.vector(suppressWarnings(MASS::fitdistr(x.resample/scale.res, densfun="gamma", start=list(shape=shape.res, scale=scale.res/scale.res), control=list(maxit=25000))$estimate))
      param.res[2] <- param.res[2] * scale.res	# Re-scale back
      D.resample[i] <- as.vector(ks.test(x.resample, "pgamma", shape=param.res[1], scale=param.res[2])$statistic)
    }
  }
  if(cdf=="pweibull") {
    if (any(x <= 0))
      stop("'x' values must be > 0 for Weibull distributions.")
    shape.x <- 1.2 / sd(log(x))
    scale.x <- exp(mean(log(x)) + 0.572/shape.x)
    # Re-scale in case the scale is much larger than 1:
    param <- as.vector(suppressWarnings(fitdistr(x/scale.x, densfun="weibull", start=list(shape=shape.x, scale=scale.x/scale.x), control=list(maxit=25000))$estimate))
    param[2] <- param[2] * scale.x	# Re-scale back
    D.obs <- as.vector(ks.test(x, "pweibull", shape=param[1], scale=param[2])$statistic)
    for (i in 1:nreps) {
      x.resample <- rweibull(n, shape=param[1], scale=param[2])
      shape.res <- 1.2 / sd(log(x.resample))
      scale.res <- exp(mean(log(x.resample)) + 0.572/shape.res)
      param.res <- as.vector(suppressWarnings(fitdistr(x.resample/scale.res, densfun="weibull", start=list(shape=shape.res, scale=scale.res/scale.res), control=list(maxit=25000))$estimate))
      param.res[2] <- param.res[2] * scale.res	# Re-scale back
      D.resample[i] <- as.vector(ks.test(x.resample, "pweibull", shape=param.res[1], scale=param.res[2])$statistic)
    }
  }
  p.value <- (sum(D.resample > D.obs) + 1) / (nreps + 1)
  se <- sqrt(p.value * (1 - p.value) / nreps)
  out <- list(D.obs=D.obs, D.resample=D.resample, p.value=p.value, se=se)
  return(out)
  }
