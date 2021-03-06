#' Lilliefors-corrected Kolmogorov-Smirnov Goodness-Of-Fit Test
#'
#' Implements the Lilliefors-corrected Kolmogorov-Smirnov test for use in
#' goodness-of-fit tests, suitable when population parameters are unknown and
#' must be estimated by sample statistics. It uses Monte Carlo simulation to
#' estimate \emph{p}-values. Using a modification of
#' \code{\link[stats]{ks.test}}, it can be used with a variety of continuous
#' distributions, including normal, lognormal, univariate mixtures of normals,
#' uniform, loguniform, exponential, gamma, and Weibull distributions. The Monte
#' Carlo algorithm can run 'in parallel.'
#'
#' @param x A numeric vector of data values (observed sample).
#' @param cdf Character string naming a cumulative distribution function. Case
#'   insensitive. Only continuous CDFs are valid. Allowed CDFs
#'   include:\itemize{\item \code{"pnorm"} for normal, \item \code{"pmixnorm"}
#'   for (univariate) normal mixture, \item \code{"plnorm"} for lognormal
#'   (log-normal, log normal), \item \code{"punif"} for uniform, \item
#'   \code{"plunif"} for loguniform (log-uniform, log uniform), \item
#'   \code{"pexp"} for exponential, \item \code{"pgamma"} for gamma, \item
#'   \code{"pweibull"} for Weibull.}
#' @param nreps Number of replicates to use in simulation algorithm.
#'   \code{Default = 4999} replicates. See \code{details} below. Should be a
#'   positive integer.
#' @param G Numeric vector of mixture components to consider, for mixture models
#'   only. \code{Default = 1:9} fits up to 9 components. Must contain positive
#'   integers. See \code{details} below.
#' @param varModel For mixture models, character string determining whether to
#'   allow equal-variance mixture components (\code{E}), variable-variance
#'   mixture components (\code{V}) or both (the \code{default}).
#' @param parallel Logical value that switches between running Monte Carlo
#'   algorithm in parallel (if \code{TRUE}) or not (if \code{FALSE}, the
#'   \code{default}).
#' @param cores Numeric value to control how many cores to build when running in
#'   parallel. \code{Default = \link[parallel]{detectCores} - 1}.
#'
#' @details The function builds a simulation distribution \code{D.sim} of length
#'   \code{nreps} by drawing random samples from the specified continuous
#'   distribution function \code{cdf} with parameters calculated from the
#'   provided sample \code{x}. Observed statistic \emph{\code{D}} and simulated
#'   test statistics are calculated using a simplified version of
#'   \code{\link[stats]{ks.test}}.
#'
#'   The default \code{nreps = 4999} provides accurate \emph{p}-values.
#'   \code{nreps = 1999} is sufficient for most cases, and computationally
#'   faster when dealing with more complicated distributions (such as univariate
#'   normal mixtures, gamma, and Weibull). See below for potentially faster
#'   parallel implementations.
#'
#'   The \emph{p}-value is calculated as the number of Monte Carlo samples with
#'   test statistics \emph{D} as extreme as or more extreme than that in the
#'   observed sample \code{D.obs}, divided by the \code{nreps} number of Monte
#'   Carlo samples. A value of 1 is added to both the numerator and denominator
#'   to allow the observed sample to be represented within the null distribution
#'   (Manly 2004); this has the benefit of avoiding nonsensical \code{p.value =
#'   0.000} and accounts for the fact that the \emph{p}-value is an estimate.
#'
#'   Parameter estimates are calculated for the specified continuous
#'   distribution, using maximum-likelihood estimates. When testing against the
#'   gamma and Weibull distributions, \code{MASS::\link[MASS]{fitdistr}} is used
#'   to calculate parameter estimates using maximum likelihood optimization,
#'   with sensible starting values. Because this incorporates an optimization
#'   routine, the simulation algorithm can be slow if using large \code{nreps}
#'   or problematic samples. Warnings often occur during these optimizations,
#'   caused by difficulties estimating sample statistic standard errors. Because
#'   such SEs are not used in the Lilliefors-corrected simulation algorithm,
#'   warnings are suppressed during these optimizations.
#'
#'   Sample statistics for the (univariate) normal mixture distribution
#'   \code{\link{pmixnorm}} are calculated using package \code{mclust}, which
#'   uses BIC to identify the optimal mixture model for the sample, and the EM
#'   algorithm to calculate parameter estimates for this model. The number of
#'   mixture components \code{G} (with default allowing up to 9 components),
#'   variance model (whether equal \code{E} or variable \code{V} variance), and
#'   component statistics (\code{mean}s, \code{sd}s, and mixing proportions
#'   \code{pro}) are estimated from the sample when calculating \code{D.obs} and
#'   passed internally when creating random Monte Carlo samples. It is possible
#'   that some of these samples may differ in their optimal \code{G} (for
#'   example a two-component input sample might yield a three-component random
#'   sample within the simulation distribution). This can be constrained by
#'   specifying that simulation BIC-optimizations only consider \code{G} mixture
#'   components.
#'
#'   Be aware that constraining \code{G} changes the null hypothesis. The
#'   default (\code{G = 1:9}) null hypothesis is that a sample was drawn from
#'   \emph{any \code{G = 1:9}-component mixture distribution}. Specifying a
#'   particular value, such as \code{G = 2}, restricts the null hypothesis to
#'   particular mixture distributions with just \code{G} components, even if
#'   simulated samples might better be represented as different mixture models.
#'
#'   The \code{LcKS(cdf = "pmixnorm")} test implements two control loops to
#'   avoid errors caused by this constraint and when working with problematic
#'   samples. The first loop occurs during model-selection for the observed
#'   sample \code{x}, and allows for estimation of parameters for the
#'   second-best model when those for the optimal model are not able to be
#'   calculated by the EM algorithm. A second loop occurs during the simulation
#'   algorithm, rejecting samples that cannot be fit by the mixture model
#'   specified by the observed sample \code{x}. Such problematic cases are most
#'   common when the observed or simulated samples have a component(s) with very
#'   small variance (i.e., duplicate observations) or when a Monte Carlo sample
#'   cannot be fit by the specified \code{G}.
#'
#'   Parellel computing can be implemented using \code{parallel = TRUE}, using
#'   the operating-system versatile \code{\link[doParallel]{doParallel-package}}
#'   and \code{\link[foreach]{foreach}} infrastructure, using a default
#'   \code{\link[parallel]{detectCores} - 1} number of cores. Parallel computing
#'   is generally advisable for the more complicated cumulative density
#'   functions (i.e., univariate normal mixture, gamma, Weibull), where maximum
#'   likelihood estimation is time-intensive, but is generally not advisable for
#'   density functions with quickly calculated sample statistics (i.e., other
#'   distribution functions). Warnings within the function provide sensible
#'   recommendations, but users are encouraged to experiment to discover their
#'   fastest implementation for their individual cases.
#'
#' @return A list containing the following components:
#'
#'   \item{D.obs}{The value of the test statistic \emph{D} for the observed
#'   sample.} \item{D.sim}{Simulation distribution of test statistics, with
#'   \code{length = nreps}. This can be used to calculate critical values; see
#'   examples.} \item{p.value}{\emph{p}-value of the test, calculated as
#'   \eqn{(\sum(D.sim > D.obs) + 1) / (nreps + 1)}.}
#'
#' @note The Kolmogorov-Smirnov (such as \code{ks.test}) is only valid as a
#'   goodness-of-fit test when the population parameters are known. This is
#'   typically not the case in practice. This invalidation occurs because
#'   estimating the parameters changes the null distribution of the test
#'   statistic; i.e., using the sample to estimate population parameters brings
#'   the Kolmogorov-Smirnov test statistic \emph{D} closer to the null
#'   distribution than it would be under the hypothesis where the population
#'   parameters are known. In other words, it is biased and results in increased
#'   Type II error rates. Lilliefors (1967, 1969) provided a solution, using
#'   Monte Carlo simulation to approximate the shape of the null distribution
#'   when the sample statistics are used to estimate population parameters, and
#'   to use this null distribution as the basis for critical values. The
#'   function \code{LcKS} generalizes this solution for a range of continuous
#'   distributions.
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}, based on
#'   code from Charles Geyer (University of Minnesota).
#'
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
#' Lc <- LcKS(x, cdf = "pnorm", nreps = 999)
#' hist(Lc$D.sim)
#' abline(v = Lc$D.obs, lty = 2)
#' print(Lc, max = 50)  # Print first 50 simulated statistics
#' # Approximate p-value (usually) << 0.05
#'
#' # Confirmation uncorrected version has increased Type II error rate when
#' #   using sample statistics to estimate parameters:
#' ks.test(x, "pnorm", mean(x), sd(x))   # p-value always larger, (usually) > 0.05
#'
#' # Confirm critical values for normal distribution are correct
#' nreps <- 9999
#' x <- rnorm(25)
#' Lc <- LcKS(x, "pnorm", nreps = nreps)
#' sim.Ds <- sort(Lc$D.sim)
#' crit <- round(c(.8, .85, .9, .95, .99) * nreps, 0)
#' # Lilliefors' (1967) critical values, using improved values from
#' #   Parsons & Wirsching (1982) (for n = 25):
#' # 0.141 0.148 0.157 0.172 0.201
#' round(sim.Ds[crit], 3)			# Approximately the same critical values
#'
#' # Confirm critical values for exponential are the same as reported by Lilliefors (1969)
#' nreps <- 9999
#' x <- rexp(25)
#' Lc <- LcKS(x, "pexp", nreps = nreps)
#' sim.Ds <- sort(Lc$D.sim)
#' crit <- round(c(.8, .85, .9, .95, .99) * nreps, 0)
#' # Lilliefors' (1969) critical values (for n = 25):
#' # 0.170 0.180 0.191 0.210 0.247
#' round(sim.Ds[crit], 3)			# Approximately the same critical values
#'
#' \dontrun{
#' # Gamma and Weibull tests require functions from the 'MASS' package
#' # Takes time for maximum likelihood optimization of statistics
#' require(MASS)
#' x <- runif(100, min = 1, max = 100)
#' Lc <- LcKS(x, cdf = "pgamma", nreps = 499)
#' Lc$p.value
#'
#' # Confirm critical values for Weibull the same as reported by Parsons & Wirsching (1982)
#' nreps <- 9999
#' x <- rweibull(25, shape = 1, scale = 1)
#' Lc <- LcKS(x, "pweibull", nreps = nreps)
#' sim.Ds <- sort(Lc$D.sim)
#' crit <- round(c(.8, .85, .9, .95, .99) * nreps, 0)
#' # Parsons & Wirsching (1982) critical values (for n = 25):
#' # 0.141 0.148 0.157 0.172 0.201
#' round(sim.Ds[crit], 3)			# Approximately the same critical values
#'
#' # Mixture test requires functions from the 'mclust' package
#' # Takes time to identify model parameters
#' require(mclust)
#' x <- rmixnorm(200, mean = c(10, 20), sd = 2, pro = c(1,3))
#' Lc <- LcKS(x, cdf = "pmixnorm", nreps = 499, G = 1:9)   # Default G (1:9) takes long time
#' Lc$p.value
#' G <- Mclust(x)$parameters$variance$G              # Optimal model has only two components
#' Lc <- LcKS(x, cdf = "pmixnorm", nreps = 499, G = G)     # Restricting to likely G saves time
#' # But note changes null hypothesis: now testing against just two-component mixture
#' Lc$p.value
#'
#' # Running 'in parallel'
#' require(doParallel)
#' set.seed(3124)
#' x <- rmixnorm(300, mean = c(110, 190, 200), sd = c(3, 15, .1), pro = c(1, 3, 1))
#' system.time(LcKS(x, "pgamma"))
#' system.time(LcKS(x, "pgamma", parallel = TRUE)) # Should be faster
#' }
#'
#' @export
#' @import mclust
#' @import MASS
#' @import doParallel
#' @import foreach
#' @import parallel
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom stats cor
#' @importFrom stats pnorm
#' @importFrom stats rnorm
#' @importFrom stats rlnorm
#' @importFrom stats runif
#' @importFrom stats rweibull
#' @importFrom stats rexp
#' @importFrom stats rgamma
#' @importFrom iterators icount
LcKS <- function(x, cdf, nreps = 4999, G = 1:9, varModel = c("E", "V"), parallel = FALSE, cores = NULL) {
  # Basic troubleshooting:
  if (missing(x) || length(x) == 0L || mode(x) != "numeric")
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x)))
    stop("'x' contains missing or infinite values")
  if (nreps < 1)
    stop("'nreps' must be a positive integer.")
  if (missing(cdf) || !is.character(cdf))
    stop("'cdf' must be a character string.")
  cdf <- tolower(cdf)
  allowed <- c("pnorm", "plnorm", "punif", "plunif", "pmixnorm", "pexp",
               "pgamma", "pweibull")
  if (!(cdf %in% allowed))
    stop("'cdf' is not a supported distribution function.")
  if (cdf == "plnorm" | cdf == "plunif" | cdf == "pexp" |
      cdf == "pgamma" | cdf == "pweibull") {
    if (any(x <= 0))
      stop("'x' values must be > 0 for your 'cdf' distribution.")
  }
  if (length(G) == 1L && G == 1 && cdf == "pmixnorm")
    stop("'G' supplied not consistent with mixture model. Use 'pnorm' or 'plnorm' instead.")
  if (!identical(G, 1:9) && cdf != "pmixnorm")
    warning("'G' is ignored except when cdf='pmixnorm'.")
  if (any(G < 1))
    stop("'G' must be a positive integer or vector of positive integers.")
  if (!parallel & !is.null(cores))
    stop("You specified to use multiple 'cores' but not to run in parallel. Set parallel = TRUE to run in parallel.")
  if (parallel & !is.null(cores)) {
    if (is.numeric(cores)) {
      if (cores < 2L | cores < 0)
        stop("Set 'cores' to integer greater than 1 to run in parallel.")
      if (!abs(cores - round(cores)) < .Machine$double.eps ^ 0.5)
        stop("Set 'cores' to integer greater than 1 to run in parallel.")
      if (cores > parallel::detectCores())
        warning("You are attempting to run this function on more cores than your computer contains. Consider reducing 'cores' to improve efficiency.")
    } else {
      stop("Set 'cores' to integer greater than 1 to run in parallel.")
    }
  }
  if (!parallel & cdf %in% c("pgamma", "pweibull", "pmixnorm"))
    warning("Consider setting 'parallel = TRUE' to increase the speed of this function.")
  if (parallel & !cdf %in% c("pgamma", "pweibull", "pmixnorm"))
    warning("Consider setting 'parallel = FALSE' to increase the speed of this function.")
  # Set up parallel infrastructure:
  if (parallel) {
    if (is.null(cores))
      cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    if (foreach::getDoParRegistered())
      message(sprintf('LcKS is running in parallel with %d worker(s) using %s [%s]\n',
                      foreach::getDoParWorkers(), foreach::getDoParName(),
                      foreach::getDoParVersion()))
  }
  # Run Monte Carlo algorithm:
  n <- length(x)
  D.sim <- rep(NA, nreps)
  if (cdf == "pnorm") {
    mean.x <- mean(x)
    sd.x <- sd(x)
    D.obs <- ks_test_stat(x, "pnorm", mean = mean.x, sd = sd.x)
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect")) %dopar% {
        x.sim <- rnorm(n, mean = mean.x, sd = sd.x)
        ks_test_stat(x.sim, "pnorm", mean = mean(x.sim), sd = sd(x.sim))
      }
    } else {
      for (i in 1:nreps) {
        x.sim <- rnorm(n, mean = mean.x, sd = sd.x)
        D.sim[i] <-
          ks_test_stat(x.sim, "pnorm", mean = mean(x.sim), sd = sd(x.sim))
      }
    }
  }
  if (cdf == "plnorm") {
    meanlog.x <- mean(log(x))
    sdlog.x <- sd(log(x))
    D.obs <- ks_test_stat(x, "plnorm", meanlog = meanlog.x, sdlog = sdlog.x)
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect")) %dopar% {
        x.sim <- rlnorm(n, meanlog = meanlog.x, sdlog = sdlog.x)
        ks_test_stat(x.sim, "plnorm", meanlog = mean(log(x.sim)),
                     sdlog = sd(log(x.sim)))
      }
    } else {
      for (i in 1:nreps) {
        x.sim <- rlnorm(n, meanlog = meanlog.x, sdlog = sdlog.x)
        D.sim[i] <- ks_test_stat(x.sim, "plnorm", meanlog = mean(log(x.sim)),
                                 sdlog = sd(log(x.sim)))
      }
    }
  }
  if (cdf == "punif") {
    min.x <- min(x)
    max.x <- max(x)
    D.obs <- ks_test_stat(x, "punif", min = min.x, max = max.x)
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect")) %dopar% {
        x.sim <- runif(n, min = min.x, max = max.x)
        ks_test_stat(x.sim, "punif", min = min(x.sim), max = max(x.sim))
      }
    } else {
      for (i in 1:nreps) {
        x.sim <- runif(n, min = min.x, max = max.x)
        D.sim[i] <-
          ks_test_stat(x.sim, "punif", min = min(x.sim), max = max(x.sim))
      }
    }
  }
  if (cdf == "plunif") {
    min.x <- min(x)
    max.x <- max(x)
    D.obs <- ks_test_stat(x, "plunif", min = min.x, max = max.x)
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect")) %dopar% {
        x.sim <- rlunif(n, min = min.x, max = max.x)
        ks_test_stat(x.sim, "plunif", min = min(x.sim), max = max(x.sim))
      }
    } else {
      for (i in 1:nreps) {
        x.sim <- rlunif(n, min = min.x, max = max.x)
        D.sim[i] <-
          ks_test_stat(x.sim, "plunif", min = min(x.sim), max = max(x.sim))
      }
    }
  }
  if (cdf == "pmixnorm") {
    G <- as.vector(G, mode = "numeric")
    varModel <- toupper(varModel)
    m <- mclust::mclustBIC(x, G = G, modelNames = varModel, verbose = FALSE)
    m <- mclust::pickBIC(m, k = sum(!is.na(m)))
    listofmod <- strsplit(names(m), ",")
    for (lom in 1:length(listofmod)) {
      mixnorm <-
        mclust::Mclust(x, G = listofmod[[lom]][2],
                       modelNames = listofmod[[lom]][1],
                       control = emControl(eps = 1e-320), verbose = FALSE)
      if (!all(is.na(mixnorm$parameters$pro)))
        break()
    }
    parameters <-
      mclust::Mclust(x, G = G, modelNames = varModel, verbose = FALSE)$parameters
    modelName <- parameters$variance$modelName
    mean.x <- parameters$mean
    sd.x <- sqrt(parameters$variance$sigmasq)
    pro.x <- parameters$pro
    D.obs <- ks_test_stat(x, "pmixnorm", mean = mean.x, sd = sd.x, pro = pro.x)
    if (parameters$variance$G == 1)
      warning("Optimal mixture model for supplied sample has a single component: it is not a mixture model. Use 'pnorm' or 'plnorm' instead.")
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect", "mclust")) %dopar% {
        repeat {
          x.sim <- rmixnorm(n, mean = mean.x, pro = pro.x, sd = sd.x)
          attempt <- try(mclust::Mclust(x.sim, G = G, modelNames = varModel,
                                        verbose = FALSE)$parameters, silent = TRUE)
          if (!inherits(attempt, "try-error")) { break }
        }
        param.sim <- attempt
        mean.sim <- param.sim$mean
        sd.sim <- sqrt(param.sim$variance$sigmasq)
        pro.sim <- param.sim$pro
        ks_test_stat(x.sim, "pmixnorm", mean = mean.sim, sd = sd.sim,
                     pro = pro.sim)
      }
    } else {
      i <- 0
      while (i < nreps) {
        x.sim <- rmixnorm(n, mean = mean.x, pro = pro.x, sd = sd.x)
        attempt <- try(mclust::Mclust(x.sim, G = G, verbose = FALSE)$parameters,
                       silent = TRUE)
        if (inherits(attempt, "try-error")) {
          next
        }
        param.sim <- attempt
        i <- i + 1
        mean.sim <- param.sim$mean
        sd.sim <- sqrt(param.sim$variance$sigmasq)
        pro.sim <- param.sim$pro
        D.sim[i] <- ks_test_stat(x.sim, "pmixnorm", mean = mean.sim,
                                 sd = sd.sim, pro = pro.sim)
      }
    }
  }
  if (cdf == "pexp") {
    rate.x <- 1 / mean(x)
    D.obs <- ks_test_stat(x, "pexp", rate = rate.x)
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect")) %dopar% {
        x.sim <- rexp(n, rate = rate.x)
        ks_test_stat(x.sim, "pexp", rate = (1 / mean(x.sim)))
      }
    } else {
      for (i in 1:nreps) {
        x.sim <- rexp(n, rate = rate.x)
        D.sim[i] <- ks_test_stat(x.sim, "pexp", rate = (1 / mean(x.sim)))
      }
    }
  }
  if (cdf == "pgamma") {
    shape.x <- mean(x) ^ 2 / var(x)
    scale.x <- var(x) / mean(x)
    # Re-scale in case the scale is much larger than 1:
    param <- as.vector(suppressWarnings(
      MASS::fitdistr(x / scale.x, densfun = "gamma",
                     start = list(shape = shape.x, scale = 1),
                     control = list(maxit = 25000))$estimate))
    param[2] <- param[2] * scale.x	# Re-scale back
    D.obs <-
      ks_test_stat(x, "pgamma", shape = param[1], scale = param[2])
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect", "MASS")) %dopar% {
        x.sim <- rgamma(n, shape = param[1], scale = param[2])
        shape.sim <- mean(x.sim) ^ 2 / var(x.sim)
        scale.sim <- var(x.sim) / mean(x.sim)
        param.sim <- as.vector(suppressWarnings(
          MASS::fitdistr(x.sim / scale.sim, densfun = "gamma",
                         start = list(shape = shape.sim, scale = 1),
                         control = list(maxit = 25000))$estimate))
        param.sim[2] <- param.sim[2] * scale.sim	# Re-scale back
        ks_test_stat(x.sim, "pgamma", shape = param.sim[1], scale = param.sim[2])
      }
    } else {
      for (i in 1:nreps) {
        x.sim <- rgamma(n, shape = param[1], scale = param[2])
        shape.sim <- mean(x.sim) ^ 2 / var(x.sim)
        scale.sim <- var(x.sim) / mean(x.sim)
        param.sim <- as.vector(suppressWarnings(
          MASS::fitdistr(x.sim / scale.sim, densfun = "gamma",
                         start = list(shape = shape.sim, scale = 1),
                         control = list(maxit = 25000))$estimate))
        param.sim[2] <- param.sim[2] * scale.sim	# Re-scale back
        D.sim[i] <-
          ks_test_stat(x.sim, "pgamma", shape = param.sim[1], scale = param.sim[2])
      }
    }
  }
  if (cdf == "pweibull") {
    shape.x <- 1.2 / sd(log(x))
    scale.x <- exp(mean(log(x)) + 0.572 / shape.x)
    # Re-scale in case the scale is much larger than 1:
    param <- as.vector(suppressWarnings(
      MASS::fitdistr(x / scale.x, densfun = "weibull",
                     start = list(shape = shape.x, scale = 1),
                     control = list(maxit = 25000))$estimate))
    param[2] <- param[2] * scale.x	# Re-scale back
    D.obs <-
      ks_test_stat(x, "pweibull", shape = param[1], scale = param[2])
    if (parallel) {
      D.sim <- foreach::foreach(iterators::icount(nreps), .combine = c,
                                .inorder = FALSE,
                                .packages = c("KScorrect", "MASS")) %dopar% {
        x.sim <- rweibull(n, shape = param[1], scale = param[2])
        shape.sim <- 1.2 / sd(log(x.sim))
        scale.sim <- exp(mean(log(x.sim)) + 0.572 / shape.sim)
        param.sim <- as.vector(suppressWarnings(
          MASS::fitdistr(x.sim / scale.sim, densfun = "weibull",
                         start = list(shape = shape.sim, scale = 1),
                         control = list(maxit = 25000))$estimate))
        param.sim[2] <- param.sim[2] * scale.sim	# Re-scale back
        ks_test_stat(x.sim, "pweibull", shape = param.sim[1], scale = param.sim[2])
        }
    } else {
      for (i in 1:nreps) {
        x.sim <- rweibull(n, shape = param[1], scale = param[2])
        shape.sim <- 1.2 / sd(log(x.sim))
        scale.sim <- exp(mean(log(x.sim)) + 0.572 / shape.sim)
        param.sim <- as.vector(suppressWarnings(
          MASS::fitdistr(x.sim / scale.sim, densfun = "weibull",
                   start = list(shape = shape.sim, scale = 1),
                   control = list(maxit = 25000))$estimate))
        param.sim[2] <- param.sim[2] * scale.sim	# Re-scale back
        D.sim[i] <-
          ks_test_stat(x.sim, "pweibull", shape = param.sim[1], scale = param.sim[2])
      }
    }
  }
  p.value <- (sum(D.sim > D.obs) + 1) / (nreps + 1)
  out <- list(D.obs = D.obs, D.sim = D.sim, p.value = p.value)
  if (parallel) parallel::stopCluster(cl)
  return(out)
}
