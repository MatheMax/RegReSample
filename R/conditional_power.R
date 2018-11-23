#' Expected conditional power
#'
#' Compute the expected conditional power weighted over all
#' minimal clinically relevant effect sizes or for a
#' point alternative
#'
#' @param n1 first stage sample size
#' @param cf stopping for futility boundary
#' @param ce stopping for efficacy boundary
#' @param n2 vector with n2-values
#' @param c2 vector with c2-values#'
#' @param weighted_alternative Should a weighted alternative be used?
#' @param delta.mcr minimal clinically relevant effect size
#'
ecp <- function(n1, cf, ce, n2, c2, weighted.alternative, delta.mcr, delta.alt,
                d.1, f.z, F.z, pi.1){
  z <- seq(cf, ce, length.out = N) # Compute nodes
  w <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes

  if(weighted.alternative == F){
    cond.pow <- function(x) { # x = c(n2, c2)
      F.z(x[2], delta.alt, x[1])
    }
    copo <- apply(cbind(n2, c2), 1, cond.pow)
  } else {
    cp <- function(x){
      cp.point <- function(delta){
        F.z(x[2], delta, x[1]) * pi.1(delta, x[3], n1)
      }
      return(h.2 * (3 * delta.alt - delta.mcr) * sum(wei.2 * sapply(d.1, cp.point)))
    }
    copo <- apply(cbind(n2, c2, z), 1, cp)
  }
  p <- h * sum(wei * w * copo)
  return(1 - p)
}



#' L_1 - norm of the derivative of n2
#'
#'
#' Approximates the L_1-norm of the derivative of n2 w.r.t. z1.
dcp.l1 <- function(n1, c2, n2, cf, ce,
                   weighted.alternative,
                   delta.mcr,
                   delta.alt,
                   d.1,
                   f.z,
                   F.z,
                   pi.1) {
  p <- cond.pow.rest(n1, cf, ce, n2, c2, weighted.alternative, delta.mcr,
                     delta.alt, d.1, F.z, pi.1)
  dis <- (ce - cf) / N
  p <- sum(abs(diff(p))) / (2 * dis)
  return(p)
}



#' Conditional Power restriction

cond.pow.rest <- function(n1, cf, ce, n2, c2, weighted.alternative, delta.mcr,
                          delta.alt, d.1, F.z, pi.1){
  z <- seq(cf, ce, length.out = N) # Compute nodes
  if(weighted.alternative == F){
    cond.pow <- function(x) { # x = c(n2, c2)
      F.z(x[2], delta.alt, x[1])
    }
    p <-  apply(cbind(n2, c2), 1, cond.pow)
  } else {
    cp <- function(x){
      cp.point <- function(delta){
        F.z(x[2], delta, x[1]) * pi.1(delta, x[3], n1)
      }
      return(h.2 * (3 * delta.alt - delta.mcr) * sum(wei.2 * sapply(d.1, cp.point)))
    }
    p <- apply(cbind(n2, c2, z), 1, cp)
  }

  return(1 - p)
}
