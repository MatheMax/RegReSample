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


ess <- function(n1, cf, ce, n2, weighted.alternative, delta.alt,
                d.2, f.z, F.z, pi.0) {
  z <- seq(cf, ce, length.out = N) # Compute nodes
  if(weighted.alternative == F){
    w <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes
    p <- h * (ce - cf) * sum(wei * w * n2)
  } else{
    ss <- function(delta){
      w <- sapply(z, function(z1) f.z(z1, delta, n1)) # Compute density at nodes
      int <- w * n2 # Multiply density by n2
      q <- h * (ce - cf) * sum(wei * w * n2) # Compute integral
      q <- q * pi.0(delta) # Multiply by prior density
    }

    p <- h.2 * (4 * delta.alt) * sum(wei.2 * sapply(d.2, ss))
  }
  p <- p + n1
  return(p)
}


#' L_1 - norm of the derivative of n2
#'
#'
#' Approximates the L_1-norm of the derivative of n2 w.r.t. z1.
#'
#' @param cf Stopping for futility boundary
#' @param ce Stopping for efficacy boundary
#' @param n2 Vector of n2-values
#'
dn2.l1 <- function(cf, ce, n2){
  #dis <- (ce - cf) / N
  p <- mean(abs(diff(n2))) #/ (2 * dis)
  return(p)
}

