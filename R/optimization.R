#' Compute the score
#'
#' score returns the weighted score
#'
#' @param n1 First-stage sample size
#' @param cf Stopping for futility boundary
#' @param ce Stopping for efficacy boundary
#' @param n2 Vector with n2-values
#' @param c2 Vector with c2-values
#' @param delta.mcr Minimal clinically relevant effect size
#' @param lambda Two-dimensional vector with penalties for CP and ||n_2^'||_1
#' @param weighted.alternative Should a weighted alternative be used?
#' @param delta.alt Point alternative effect size if weighted.alternative = F,
#' prior mean otherwise
#' @param tau Standard deviation of prior density
#'
score <- function(n1, cf, ce, n2, c2, lambda,
                  weighted.alternative, delta.mcr, delta.alt, tau) {
  p <- ess(n1, cf, ce, n2, weighted.alternative, delta.alt, tau)
  if(lambda[1] != 0){
    p <- p - lambda[1] *
      ecp(n1, cf, ce, n2, c2, delta.mcr, weighted.alternative, delta.alt, tau)
  }
  if(lambda[2] != 0){
    p <- p  + lambda[2] * dn2.l1(cf, ce, n2)
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
#' @param lambda Two-dimensional vector with penalties for CP and ||n_2^'||_1
#' @param delta.mcr Minimal clinically relevant effect size
#' @param weighted.alternative Should a weighted alternative be used?
#' @param delta.alt Point alternative effect size if weighted.alternative = F,
#' prior mean otherwise
#' @param tau Standard deviation of the prior density
#' @param n.max Maximal sample size per stage
#'
#' Smoothing splines are used to approximate n_2 and c_2 functions.
#' All integrals are approximates by Boole's Newton-Cotes formula
#' or by the \code{R}-routine \code{integrate}.
#'
#' @export

opt_design <- function(alpha,
                beta,
                lambda,
                weighted.alternative = FALSE,
                delta.mcr = 0,
                delta.alt = .3,
                tau = .1,
                n.max = Inf) {

  start <- c(50, 0, 2, seq(50, 10, length.out = N), seq(2, 0, length.out = N))
  low <- c(0, -1, qnorm(1 - alpha), rep(0, N), rep(-1, N))
  up <- c(n.max, 4, 4, rep(n.max, N), rep(4, N))

  optimum <- nloptr::nloptr(
    x0          = start,
    eval_f      = function(x) {score(x[1], x[2], x[3],
                                                 x[4 : (N+3)],
                                                 x[(N+4) : length(start)],
                                                 lambda,
                                                 weighted.alternative,
                                                 delta.mcr,
                                                 delta.alt,
                                                 tau)
                },
    eval_g_ineq = function(x) { return( c(
                  x[2] - x[3] + 0.1,
                  toe(x[2], x[3], x[(N+4) : length(start)], x[1]) - alpha,
                  1 - beta - pow(x[1], x[2], x[3],
                                 x[4 : (N+3)],
                                 x[(N+4) : length(start)],
                                 weighted.alternative,
                                 delta.mcr,
                                 delta.alt,
                                 tau)
    ) )
    },
    lb = low,
    ub = up,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.0001,
      maxeval = 99999,
      maxtime = 16200
    )
  )

  n1 <- optimum$solution[1]
  cf <- optimum$solution[2]
  ce <- optimum$solution[3]
  n2 <- optimum$solution[4:(N+3)]
  c2 <- optimum$solution[(N+4):length(start)]
  x <- seq(cf, ce, length.out = N)


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
