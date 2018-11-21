#' L_1 - norm of the derivative of n2
#'
#' Approximates the L_1-norm of the derivative of n2 w.r.t. z1.
#'
#' @param cf Stopping for futility boundary
#' @param ce Stopping for efficacy boundary
#' @param n2 Vector of n2-values
#'
dn2.l1 <- function(cf, ce, n2){
  dis <- (ce - cf) / length(n2)
  p <- sum(abs(diff(n2))) / (2 * dis)
  return(p)
}

# Conditional power
dcp.l1 <- function(c2, n2, cf, ce,
                   weighted.alternative = F,
                   delta.mcr = 0,
                   delta.alt = .3,
                   tau = .1) {
  dis <- (ce - cf) / length(n2)
  z <- seq(cf, ce, length.out = 10)
  if(weighted.alternative == T){
    p <- sapply(z, function(z1) bcp(z1, c2, n2, cf, ce, delta.mcr, delta.alt, tau))
  } else{
    p <- sapply(z, function(z1) cp(z1, c2, n2, cf, ce, delta.alt))
  }
  p <- sum(abs(diff(p))) / (2 * dis)
  return(p)
}
