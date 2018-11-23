#' Compute the score
#'
#' score returns the weighted score
#'
#' @param n1 First-stage sample size
#' @param cf Stopping for futility boundary
#' @param ce Stopping for efficacy boundary
#' @param n2 Vector with n2-values
#' @param c2 Vector with c2-values
#' @param lambda Three-dimensional vector with penalties for ECP, ||n_2^'||_1, ||CP^'||_1
#' @param weighted.alternative Should a weighted alternative be used?
#' @param delta.mcr Minimal clinically relevant effect size
#' @param delta.alt Point alternative effect size if weighted.alternative = F,
#' prior mean otherwise
#' @param tau Standard deviation of prior density
#'
score <- function(n1, cf, ce, n2, c2, lambda,
                  weighted.alternative, delta.mcr, delta.alt,
                  d.1, d.2, f.z, F.z, pi.0, pi.1) {
  p <- ess(n1, cf, ce, n2, weighted.alternative, delta.alt, d.2, f.z, F.z, pi.0)
  if(lambda[1] != 0){
    p <- p - lambda[1] *
      ecp(n1, cf, ce, n2, c2, weighted.alternative, delta.mcr, delta.alt,
          d.1, f.z, F.z, pi.1)
  }
  if(lambda[2] != 0){
    p <- p  + lambda[2] * dn2.l1(cf, ce, n2)
  }
  if(lambda[3] != 0){
    p <- p  + lambda[3] * dcp.l1(n1, c2, n2, cf, ce, weighted.alternative,
                                 delta.mcr, delta.alt, d.1, f.z, F.z, pi.1)
  }

  return(p)
}


#' Compute an optimal design
#'
#' \code{opt_design} returns a design that is optimal for
#' specified parameters and a desired mix of expected sample size,
#' its variability and conditional power and its variability
#'
#' @param alpha Maximal type one error rate
#' @param beta Maximal type two error rate
#' @param beta.2 Conditional power should be at least 1 - beta.2
#' @param lambda Three-dimensional vector with penalties for ECP, ||n_2^'||_1, ||CP^'||_1
#' @param prior Should a normal or a point prior be used?
#' @param power Should the prior density or predictive power be used for conditional power?
#' @param delta.mcr Minimal clinically relevant effect size
#' @param delta.alt Point alternative effect size if weighted.alternative = F,
#' prior mean otherwise
#' @param tau Standard deviation of the prior density
#' @param n.max Maximal sample size per stage
#' @param distribution Is the test statistic normal or t-distributed?
#'
#' Smoothing splines are used to approximate n_2 and c_2 functions.
#' All integrals are approximates by Boole's Newton-Cotes formula.
#'
#' @export

opt_design <- function(alpha,
                beta,
                beta.2,
                lambda,
                prior = c("normal", "point"),
                power = c("equal", "predictive"),
                delta.mcr = 0,
                delta.alt = .3,
                tau = .1,
                n.max = Inf,
                distribution = c("normal", "t")
                ) {

  # Build sequence of effect sizes in order to compute integrals
  d.1 <- seq(delta.mcr, 3 * delta.alt, length.out = N.2)
  d.2 <- seq(-delta.alt, 3 * delta.alt, length.out = N.2)

  # Define distribution of z_1
  if(distribution == "normal"){
    f.z <- function(z, delta, n){
    dnorm(z, mean = sqrt(n) * delta, sd = 1)  # Lebesgue density of Z_1
    }
    F.z <- function(z, delta, n) {
    pnorm(z, mean = sqrt(n) * delta, sd = 1)  # Distribution function of Z_1
    }
  } else if(distribution == "t"){
    f.z <- function(z, delta, n){
    dt(z, df = n - 1, ncp = sqrt(n) * delta) # Lebesgue density of Z_1
    }
    F.z <- function(z, delta, n) {
    pt(z, df = n - 1, ncp = sqrt(n) * delta)  # Distribution function of Z_1
    }
  }

  # Define prior
  if(prior == "point"){
  weighted.alternative = F
  pi.0 <- function(delta){
    dnorm(delta, mean = delta.alt, sd = tau)
  }
  } else if (prior == "normal"){
  weighted.alternative = T
  pi.0 <- function(delta){
    dnorm(delta, mean = delta.alt, sd = tau)
  }
  }

  # Define power
  if(power == "equal"){
  pi.1 <- function(delta, z1, n1){
    dnorm(delta, mean = delta.alt, sd = tau)
  }
  } else if(power == "predictive"){
  pi.1 <- function(delta, z1, n1){
    t <- 1 / (n1 + 1 / tau^2)
    dnorm(delta, mean = t * (delta.alt / tau^2 + sqrt(n1) * z1), sd = sqrt(t))
  }
  }


  # Define start values and boundaries
  start <- c(min(50, n.max - 1), 0, min(2, qnorm(1 - alpha)),
             seq(min(50, n.max - 1), min(10, n.max - 1), length.out = N),
             seq(2, 0, length.out = N))
  low <- c(0, -1, qnorm(1 - alpha),
           rep(0, N), rep(-1, N))
  up <- c(n.max, 4, 4,
          rep(n.max, N), rep(4, N))

  # Define errors
  if(pow(up[1], low[2], up[3], up[4 : (N+3)], low[(N+4) : length(start)],
         weighted.alternative, delta.mcr, delta.alt, d.1, f.z, F.z,
         pi.0) < 1 - beta){
    return("The maximal sample size is too low for the desired power!")
    break
  }

  if(lambda[1] < 0 | lambda[2] < 0 | lambda[3] < 0){
    return("lambda must not be negative!")
    break
  }

  # Optimize
  optimum <- nloptr::nloptr(
    x0          = start,
    eval_f      = function(y) {score(y[1], y[2], y[3],
                                                 y[4 : (N+3)],
                                                 y[(N+4) : length(start)],
                                                 lambda,
                                                 weighted.alternative,
                                                 delta.mcr,
                                                 delta.alt,
                                                 d.1,
                                                 d.2,
                                                 f.z,
                                                 F.z,
                                                 pi.0,
                                                 pi.1)
                },
    eval_g_ineq = function(y) { return( c(
                  y[2] - y[3] + 0.1,
                  diff(y[(N+4) : length(start)]),
                  toe(y[1], y[2], y[3], y[(N+4) : length(start)], f.z, F.z) - alpha,
                  1 - beta - pow(y[1], y[2], y[3],
                                 y[4 : (N+3)],
                                 y[(N+4) : length(start)],
                                 weighted.alternative,
                                 delta.mcr,
                                 delta.alt,
                                 d.1,
                                 f.z,
                                 F.z,
                                 pi.0),
                  rep(1 - beta.2, N) - cond.pow.rest(y[1], y[2], y[3],
                                                   y[4 : (N+3)],
                                                   y[(N+4) : length(start)],
                                                   weighted.alternative,
                                                   delta.mcr,
                                                   delta.alt,
                                                   d.1,
                                                   F.z,
                                                   pi.1)
    ) )
    },
    lb = low,
    ub = up,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.00001,
      maxeval = 99999,
      maxtime = 16200
    )
  )

  n1 <- optimum$solution[1]
  cf <- optimum$solution[2]
  ce <- optimum$solution[3]
  n2 <- optimum$solution[4 : (N + 3)]
  c2 <- optimum$solution[(N + 4) : length(start)]
  x <- seq(cf, ce, length.out = N)

  y <- optimum$solution

  if(toe(y[1], y[2], y[3], y[(N+4) : length(start)], f.z, F.z) - alpha > .005){
    return("The chosen parameters do not allow to protect the type one error rate!")
    break
  }

  if(1 - beta - pow(y[1], y[2], y[3], y[4 : (N+3)], y[(N+4) : length(start)],
                    weighted.alternative, delta.mcr, delta.alt,
                    d.1, f.z, F.z, pi.0) > .01){
    return("The chosen parameters do not allow to achieve the desired power!")
    break
  }


  n2_out <- function(z) {
    spl <- smooth.spline(x, n2)
    p <- ifelse(cf <= z && z <= ce, predict(spl, z, 0)$y, 0)
    return(p)
  }

  c2_out <- function(z) {
    spl <- smooth.spline(x, c2)
    if(z < cf) {
      p = Inf
    }else if(z > ce) {
      p = -Inf
    } else {
      p <- predict(spl, z, 0)$y
    }
    return(p)
  }

  d <- design(cf, ce, n1, n2_out, c2_out)
  return(d)
}
