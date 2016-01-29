#' @rdname dmixnorm
#' @export
pmixnorm <- function(q, mean, sd, pro) {
  if(mode(q) != "numeric")
    stop("'q' must be a non-empty numeric vector")
  if(any(missing(mean), missing(sd)))
    stop("'mean' and 'sd' not provided, without default.")
  mean <- as.vector(mean, mode="numeric")
  G <- length(mean)
  sd <- as.vector(sd, mode="numeric")
  if (missing(pro)) {
    pro <- rep(1/G, G)
    warning("mixing proportion 'pro' not provided. Assigned equal proportions
            by default.")
  }
  if(any(pro < 0L, sd < 0L))
    stop("'pro' and 'sd' must not be negative.")
  lpro <- length(pro)
  if(length(sd)==1) sd[seq(G)] <- sd[1]
  lsd <- length(sd)
  if(G < lsd | G < lpro | (lsd > 1 & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode="numeric")
  pro <- pro/sum(pro)
  cdf <- rep(0, length(q))
  for(g in seq(G)) { cdf <- cdf + pro[g] * pnorm(q, mean[g], sd[g]) }
  return(cdf)
}
