#' Expected sample size
#'
#' Computes the expected sample size under a prespecified (average of)
#' effect sizes
#'
#' @param n1 First stage sample size
#' @param cf stoppinf for futility boundary
#' @param ce Stopping for efficacy boundary
#' @param n2 Vector with n2-values
#' @param weighted_alternative Should a weighted average be used?


ess <- function(n1, cf, ce, n2, weighted.alternative, delta.alt) {
  z <- seq(cf, ce, length.out = length(n2)) # Compute nodes
  if(weighted.alternative == FALSE){
    w <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes
    p <- h * (ce - cf) * sum(wei * w * n2)
  } else{
    p <- integrate(Vectorize(function(delta) {
      w <- sapply(z, function(z1) f.z(z1, delta, n1)) # Compute density at nodes
      int <- w * n2 # Multiply density by n2
      q <- h * (ce - cf) * sum(wei * w * n2) # Compute integral
      q <- q * pi_0(delta, delta.alt) # Multiply by prior density
    }),
    -Inf,
    Inf)$value
  }
  p <- p + n1
  return(p)
}
