fixed_design <- function(alpha, beta, delta.mcr){
  c <- qnorm(1 - alpha) # Ensure type one error rate

  pow <- function(n){
    integrate(function(delta){
      pi_0(delta) * pnorm(sqrt(n) * delta - c)
    },
    delta.mcr,
    Inf)$value
  } # Define power function

  optimum <- nloptr::nloptr(
    x0 = 100,
    eval_f      = function(x) x,
    eval_g_ineq = function(x) 1 - beta - pow(x),
    lb = 0,
    ub = Inf,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.0001,
      maxeval = 99999,
      maxtime = 16200
    )
  ) # Compute n

  n <- optimum$solution

  return(c(n, c))

}
