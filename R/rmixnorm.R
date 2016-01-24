#' @rdname dmixnorm
#' @export
#' @importFrom stats rnorm
rmixnorm <- function (n, mean, sd, pro) {
  if(mode(n) != "numeric" | n <= 1L)
    stop("'n' must be a positive, non-empty numeric vector.")
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
  if(length(sd)==1) ( modelName <- "E" ) else ( modelName <- "V" )
  if(length(sd)==1) sd[seq(G)] <- sd[1]
  lsd <- length(sd)
  if(G < lsd | G < lpro | (lsd > 1 & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode="numeric")
  pro <- pro/sum(pro)
  clabels <- sample(1:G, size=n, replace=TRUE, prob=pro)
  ctable <- tabulate(clabels, nbins=G)
  x <- rep(0, n)
  for (k in 1:G) { x[clabels == k] <- mean[k] + rnorm(ctable[k], sd=sd[k]) }
  structure(as.vector(x), modelName = modelName)
}
