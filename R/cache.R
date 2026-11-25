# Cache of simulated null distributions, used by LcKS(cache = TRUE).
#
# When 'cdf' is a location-scale family (possibly after taking logs) and its
# parameters are estimated equivariantly, the null distribution of D depends
# only on the sample size, so that one simulation can serve every sample of the
# same size. For the Weibull distribution, this holds up to the precision of
# MASS::fitdistr().

null_cache <- new.env(parent = emptyenv())

# Quantile functions of the standard distributions. The names of this list
# define the 'cdf' values that can be cached.
#' @importFrom stats qnorm
#' @importFrom stats qlnorm
#' @importFrom stats qunif
#' @importFrom stats qexp
#' @importFrom stats qweibull
cacheable_cdfs <- list(
  pnorm = qnorm,
  plnorm = qlnorm,
  punif = qunif,
  plunif = function(p) qlunif(p, min = 1, max = exp(1)),
  pexp = qexp,
  pweibull = function(p) qweibull(p, shape = 1)
)

# Simulated null distribution of D for samples of size 'n', computed by LcKS()
# on its first use and read from the cache afterwards. The simulation is seeded
# so that its result depends only on the arguments, and leaves the user's
# random number stream untouched.
#' @importFrom stats ppoints
cached_D_sim <- function(cdf, n, nreps) {
  key <- paste(cdf, n, nreps, sep = "_")
  if (is.null(null_cache[[key]])) {
    # Deterministic sample: any sample of size 'n' has the same null distribution
    x.std <- cacheable_cdfs[[cdf]](ppoints(n))
    null_cache[[key]] <- with_fixed_seed(
      90717 + n,
      suppressWarnings(LcKS(x.std, cdf = cdf, nreps = nreps)$D.sim)
    )
  }
  null_cache[[key]]
}

# Evaluate 'expr' with a fixed seed and the default generators of R, then
# restore the user's generators and seed.
with_fixed_seed <- function(seed, expr) {
  old.kind <- RNGkind()
  old.seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    # sample.kind = "Rounding" triggers a warning
    suppressWarnings(RNGkind(old.kind[1], old.kind[2], old.kind[3]))
    if (is.null(old.seed))
      rm(".Random.seed", envir = globalenv())
    else
      assign(".Random.seed", old.seed, envir = globalenv())
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}
