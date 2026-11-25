# Tests of LcKS(cache = TRUE)
library(KScorrect)

null_cache <- KScorrect:::null_cache
clear <- function() rm(list = ls(null_cache), envir = null_cache)
warns <- function(expr) tryCatch({ expr; FALSE }, warning = function(w) TRUE)
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

# Pairs of samples with very different parameters
set.seed(1)
n <- 10
samples <- list(pnorm = list(rnorm(n), rnorm(n, 50, 9)),
                plnorm = list(rlnorm(n), rlnorm(n, 3, 2)),
                punif = list(runif(n), runif(n, -40, 90)),
                plunif = list(rlunif(n, 1, 2), rlunif(n, 0.01, 900)),
                pexp = list(rexp(n), rexp(n, 0.01)))
stopifnot(setequal(names(samples), names(KScorrect:::cacheable_cdfs)))

for (cdf in names(samples)) {
  x <- samples[[cdf]][[1]]
  y <- samples[[cdf]][[2]]

  # The null distribution does not depend on the estimated parameters
  set.seed(2)
  a <- LcKS(x, cdf, nreps = 199)
  set.seed(2)
  b <- LcKS(y, cdf, nreps = 199)
  stopifnot(all.equal(a$D.sim, b$D.sim))

  # The result depends neither on the state of the cache nor on the seed
  clear()
  cold <- LcKS(x, cdf, nreps = 199, cache = TRUE)
  set.seed(3)
  warm <- LcKS(x, cdf, nreps = 199, cache = TRUE)
  stopifnot(identical(cold, warm))
  clear()
  first <- LcKS(y, cdf, nreps = 199, cache = TRUE)
  stopifnot(identical(first$D.sim, cold$D.sim),
            identical(LcKS(x, cdf, nreps = 199, cache = TRUE), cold))

  # Same statistic as without the cache, and consistent p-value
  stopifnot(identical(cold$D.obs, a$D.obs),
            length(cold$D.sim) == 199,
            cold$p.value == (sum(cold$D.sim > cold$D.obs) + 1) / 200)

  # Same null distribution as without the cache
  a <- LcKS(x, cdf, nreps = 4999)
  b <- LcKS(x, cdf, nreps = 4999, cache = TRUE)
  stopifnot(abs(quantile(a$D.sim, 0.95) - quantile(b$D.sim, 0.95)) < 0.01)
}

# One entry per distribution, sample size and number of replicates
clear()
x <- rexp(12)
invisible(LcKS(x, "pexp", nreps = 99, cache = TRUE))
invisible(LcKS(2 * x, "PEXP", nreps = 99, cache = TRUE))
stopifnot(length(ls(null_cache)) == 1)
invisible(LcKS(x[-1], "pexp", nreps = 99, cache = TRUE))
invisible(LcKS(x, "pexp", nreps = 199, cache = TRUE))
invisible(LcKS(x, "pnorm", nreps = 99, cache = TRUE))
stopifnot(length(ls(null_cache)) == 4)

# The random number stream is left untouched, whatever the generators
reference <- LcKS(x, "pnorm", nreps = 99, cache = TRUE)
for (kind in list(c("Mersenne-Twister", "Inversion", "Rejection"),
                  c("L'Ecuyer-CMRG", "Ahrens-Dieter", "Rounding"))) {
  suppressWarnings(RNGkind(kind[1], kind[2], kind[3]))
  set.seed(4)
  before <- .Random.seed
  clear()
  stopifnot(identical(LcKS(x, "pnorm", nreps = 99, cache = TRUE), reference),
            identical(.Random.seed, before),
            identical(RNGkind(), kind))
}
rm(.Random.seed)
clear()
invisible(LcKS(x, "pnorm", nreps = 99, cache = TRUE))
stopifnot(!exists(".Random.seed", envir = globalenv(), inherits = FALSE),
          identical(RNGkind(), kind))
RNGkind("default", "default", "default")
set.seed(5)

# No cache by default, nor for the other distributions
clear()
invisible(LcKS(x, "pexp", nreps = 99))
stopifnot(warns(LcKS(x, "pgamma", nreps = 1, cache = TRUE)))
g <- suppressWarnings(LcKS(x, "pgamma", nreps = 9, cache = TRUE))
stopifnot(length(g$D.sim) == 9, length(ls(null_cache)) == 0)

# Incompatible or invalid arguments
stopifnot(warns(LcKS(x, "pexp", nreps = 99, cache = TRUE, parallel = TRUE)),
          fails(LcKS(x, "pexp", nreps = 99, cache = NA)),
          fails(LcKS(x, "pexp", nreps = 99, cache = "yes")))

# Published critical values (n = 25; levels 0.8, 0.85, 0.9, 0.95, 0.99)
crit <- function(cdf) {
  D.sim <- LcKS(rexp(25), cdf, nreps = 9999, cache = TRUE)$D.sim
  unname(quantile(D.sim, c(0.8, 0.85, 0.9, 0.95, 0.99)))
}
# Lilliefors (1967), improved by Parsons and Wirsching (1982)
stopifnot(abs(crit("pnorm") - c(0.141, 0.148, 0.157, 0.172, 0.201)) < 0.005)
# Lilliefors (1969)
stopifnot(abs(crit("pexp") - c(0.170, 0.180, 0.191, 0.210, 0.247)) < 0.005)
