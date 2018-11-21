#' Type one error rate
#'
#' Compute the type one error for given decision boundaries
#'
#' @param cf stopping for futility boundary
#' @param ce stopping for efficacy boundary
#' @param c2 vector with c2-values
#'
toe <- function(cf, ce, c2, n1) {
  p <- pnorm(-cf)
  z <- seq(cf, ce, length.out = length(c2))
  int <- h *  (ce - cf) * sum(wei * F.z(c2, 0, n1) * f.z(z, 0, n1))
  #int <- h *  (ce - cf) * sum(wei * pnorm(c2) * dnorm(z))
  p <- p - int
  return(p)
}

#' Power
#'
#' Compute the power of a design
#'
#' @param n1 First stage sample size
#' @param cf stoppinf for futility boundary
#' @param ce Stopping for efficacy boundary
#' @param n2 Vector with n2-values
#' @param c2 Vector with c2-values
#' @param delta.mcr Minimal clinically relevant effect size
#' @param weighted_alternative Should a weighted alternative be used?
#'
pow <- function(n1, cf, ce, n2, c2, weighted.alternative, delta.mcr, delta.alt, tau) {
  z <- seq(cf, ce, length.out = length(c2))
  if(weighted.alternative == FALSE) {
    fs <- pnorm(sqrt(n1) * delta.alt - cf)
    den <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes
    cond.pow <- function(x) {
      F.z(x[2], delta.alt, x[1])  # x = c(n2, c2)
    }
    int <- den * apply(cbind(n2, c2), 1, cond.pow)
    int <- h * (ce - cf) * sum(wei * int)
    pow <- (fs - int)
  } else {
    pow <- integrate(Vectorize(function(delta) {
      fs <- pnorm(sqrt(n1) * delta - cf)
      den <- sapply(z, function(z1) f.z(z1, delta, n1)) # Compute density at nodes
      cond.pow <- function(x) {
        F.z(x[2], delta, x[1]) # x = c(n2, c2)
      }
      int <- den * apply(cbind(n2, c2), 1, cond.pow)
      int <- h * (ce - cf) * sum(wei * int)
      q <- (fs - int) * pi.0(delta, delta.alt, tau)
      return(q)
  }),
  delta.mcr,
  Inf)$value
  }
  return(pow)
}


#' Conditional Power restriction

cond.pow.rest <- function(n1, cf, ce, n2, c2, delta.mcr, weighted.alternative, delta.alt, tau){
  z <- seq(cf, ce, length.out = length(n2)) # Compute nodes
  if(weighted.alternative == FALSE){
    cond.pow <- function(x) { # x = c(n2, c2)
      F.z(x[2], delta.alt, x[1])
    }
    p <-  apply(cbind(n2, c2), 1, cond.pow)
  } else {
    cp <- function(x){
      integrate(Vectorize(function(delta){
        copo <-  F.z(x[2], delta, x[1]) * pi.1(delta, delta.alt, tau, x[3], n1)
        q <-  copo
        return(q)
    }),
    delta.mcr,
    Inf)$value
    }

    p <- apply(cbind(n2, c2, z), 1, cp)
  }

  return(1 - p)
}
