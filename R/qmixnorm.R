#' @rdname dmixnorm
#' @export
qmixnorm <- function (npts=100, nreps=5000, mean, sd, pro) {
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
  pps <- matrix(NA, nreps, npts)
  for(p in 1:nreps) { pps[p,] <- sort(rmixnorm(npts, mean=mean, pro=pro, sd=sd)) }
  exp <- apply(pps, 2, mean)
  return(exp)
  }
