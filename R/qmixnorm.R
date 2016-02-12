#' @rdname dmixnorm
#' @export
#' @importFrom stats quantile
qmixnorm <- function (p, mean, sd, pro, nr = 200000)  {
  if(mode(p) != "numeric")
    stop("'p' must be a non-empty numeric vector")
  if(mode(nr) != "numeric" | nr <= 1L)
    stop("'nr' must be a positive, non-empty numeric vector.")
  if (any(missing(mean), missing(sd)))
    stop("'mean' and 'sd' not provided, without default.")
  mean <- as.vector(mean, mode = "numeric")
  G <- length(mean)
  sd <- as.vector(sd, mode = "numeric")
  if (missing(pro)) {
    pro <- rep(1/G, G)
    warning("mixing proportion 'pro' not provided. Assigned equal proportions by default.")
  }
  if (any(pro < 0L, sd < 0L))
    stop("'pro' and 'sd' must not be negative.")
  lpro <- length(pro)
  if(length(sd)==1) sd[seq(G)] <- sd[1]
  lsd <- length(sd)
  if(G < lsd | G < lpro | (lsd > 1 & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode = "numeric")
  pro <- pro/sum(pro)
  samp <- rmixnorm(nr, mean = mean, sd = sd, pro = pro)
  quants <- stats::quantile(samp, p)
  if (any(p < .01, p > 0.99))
    warning("quantile approximations may not be reliable for probability values close to 0 or 1.")
  tol <- -log10(.Machine$double.eps)
  if(any(round(p, tol)==0L)) quants[which(round(p, tol)==0L)] <- -Inf
  if(any(round(p, tol)==1L)) quants[which(round(p, tol)==1L)] <- Inf
  return(as.vector(quants))
}
