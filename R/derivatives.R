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
  p <- 0
  for(i in 2:length(n2)) {
    p <- p + abs(n2[i] - n2[i-1])
  }
  p <- p / (2 * dis)
  return(p)
}

# Conditional power
