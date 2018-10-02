#' Type one error rate
#'
#' Compute the type one error for given decision boundaries
#'
#' @param cf stopping for futility boundary
#' @param ce stopping for efficacy boundary
#' @param c2 vector with c2-values
#'
toe <- function(cf, ce, c2) {
  p <- pnorm(-cf)
  z <- seq(cf, ce, length.out = length(c2))
  int <- h *  (ce - cf) * sum(wei * pnorm(c2) * dnorm(z))
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
pow <- function(n1, cf, ce, n2, c2, delta.mcr, weighted.alternative, delta.alt) {
  z <- seq(cf, ce, length.out = length(c2))
  if(weighted.alternative == FALSE) {
    fs <- pnorm(sqrt(n1) * delta.alt - cf)
    den <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes
    cond.pow <- function(x) {
      F.z(x[2], delta.alt, x[1])
      #pnorm(x[2] - sqrt(x[1]) * delta.alt) # x = c(n2, c2)
    }
    int <- den * apply(cbind(n2, c2), 1, cond.pow)
    int <- h * (ce - cf) * sum(wei * int)
    pow <- (fs - int)
  } else {
    pow <- integrate(Vectorize(function(delta) {
      fs <- pnorm(sqrt(n1) * delta - cf)
      den <- sapply(z, function(z1) f.z(z1, delta, n1)) # Compute density at nodes
      cond.pow <- function(x) {
        F.z(x[2], delta, x[1])
        #pnorm(x[2] - sqrt(x[1]) * delta) # x = c(n2, c2)
      }
      int <- den * apply(cbind(n2, c2), 1, cond.pow)
      int <- h * (ce - cf) * sum(wei * int)
      q <- (fs - int) * pi_0(delta, delta.alt)
      return(q)
  }),
  delta.mcr,
  Inf)$value
  }
  return(pow)
}
