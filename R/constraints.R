#' Type one error rate
#'
#' Compute the type one error for given decision boundaries
#'
#' @param cf stopping for futility boundary
#' @param ce stopping for efficacy boundary
#' @param c2 vector with c2-values
#'
toe <- function(n1, cf, ce, c2, f.z, F.z) {
  p <- F.z(-cf, 0, n1)
  z <- seq(cf, ce, length.out = N)
  int <- h *  (ce - cf) * sum(wei * F.z(c2, 0, n1) * f.z(z, 0, n1))
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
pow <- function(n1, cf, ce, n2, c2, weighted.alternative,
                delta.mcr, delta.alt, d.1, f.z, F.z, pi.0, pi.1) {
  z <- seq(cf, ce, length.out = N)
  if(weighted.alternative == F) {
    fs <- 1 - F.z(cf,delta.alt, n1)
    den <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes
    cond.pow <- function(x) {
      F.z(x[2], delta.alt, x[1])  # x = c(n2, c2)
    }
    int <- den * apply(cbind(n2, c2), 1, cond.pow)
    int <- h * (ce - cf) * sum(wei * int)
    pow <- (fs - int)
  } else {
    pow.point <- function(delta) {
      fs <- 1 - F.z(cf, delta, n1)
      den <- sapply(z, function(z1) f.z(z1, delta, n1)) # Compute density at nodes
      cond.pow <- function(x) {
        F.z(x[2], delta, x[1]) # x = c(n2, c2)
      }
      int <- den * apply(cbind(n2, c2), 1, cond.pow)
      int <- h * (ce - cf) * sum(wei * int)
      q <- (fs - int) * pi.0(delta)
      return(q)
    }
    pow <- h.2 * (3 * delta.alt - delta.mcr) * sum(wei.2 * sapply(d.1, pow.point))
  }

  return(pow)
}

