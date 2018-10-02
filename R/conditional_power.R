#' Conditional power for specific delta
#'
#' Compute the conditional power for a specified effect size
#' and speficied parameters
#'
#' @param z1 the observed value of the first stage test statistic
#' @param c2 a vector with c2-values
#' @param n2 a vector with n2-values
#' @param cf stopping for futility boundary
#' @param ce stopping for efficacy boundary
#' @param delta effect size
#'
cp <- function(z1, c2, n2, cf, ce, delta) {
  z <- seq(cf, ce, length.out = length(n2))
  c <- splinefun(z, c2)
  n <- splinefun(z, n2)
  p <- F(c(z1), delta, n(z1))
  return(p)
}


#' Weighted conditional power over all clinically relevant deltas
#'
#' Compute the conditional power for a weighted average of
#' all delta-values which are clinically relevant
#'
#' @param z1 the observed value of the first stage test statistic
#' @param c2 a vector with c2-values
#' @param n2 a vector with n2-values
#' @param cf stopping for futility boundary
#' @param ce stopping for efficacy boundary
#' @param delta.mcr minimal clinically relevant effect
#'
bcp <- function(z1, c2, n2, cf, ce, delta.mcr, delta.alt, tau) {
  z <- seq(cf, ce, length.out = length(n2))
  c <- splinefun(z, c2)
  n <- splinefun(z, n2)
  p <- integrate(Vectorize(function(delta) {
              F(c(z1), delta, n(z1)) * pi.1(delta, delta.alt, tau, z1, n1)
    #evtl kann man hier laufzeit gewinnen
            },
            delta.mcr,
            1))
  return(p$value)
}


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
#' @param c2 vector with c2-values
#' @param delta.mcr minimal clinically relevant effect size
#' @param weighted_alternative Should a weighted alternative be used?
#'
ecp <- function(n1, cf, ce, n2, c2, delta.mcr, weighted.alternative, delta.alt, tau){
  z <- seq(cf, ce, length.out = length(n2)) # Compute nodes
  if(weighted.alternative == FALSE){
    w <- sapply(z, function(z1) f.z(z1, delta.alt, n1)) # Compute density at nodes
    cond.pow <- function(x) { # x = c(n2, c2)
      F.z(x[2], delta.alt, x[1])
    }
    copo <- apply(cbind(n2, c2), 1, cond.pow)
    p <- h * (ce - cf) * sum(wei * w * copo)
  } else {
    p <- integrate(Vectorize(function(delta){
      w <- sapply(z, function(z1) f.z(z1, delta, n1)) # Compute density at nodes
      cond.pow <- function(x) { # x = c(n2, c2, z1)
        F.z(x[2], delta, x[1]) * pi.1(delta, delta.alt, tau, x[3], n1)
      }
      copo <- apply(cbind(n2, c2, z), 1, cond.pow)
      q <- h * (ce - cf) * sum(wei * w * copo)
      return(q)
      }),
      delta.mcr,
      Inf)$value
  }
  p <- p / (ce - cf)
  return(p)
}
