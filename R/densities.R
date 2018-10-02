# Prior density
pi.0 <- function(delta, delta.alt, tau){
  dnorm(delta, mean = delta.alt, sd = tau)
}

# Posterior density
pi.1 <- function(delta, delta.alt, tau, z1, n1){
  dnorm(delta, mean = delta.alt, sd = tau)
  #For predictive power use:
  #t <- 1 / (n1 + 1 / tau^2)
  #dnorm(delta, mean = t * (delta.alt / tau^2 + sqrt(n1) * z1), sd = sqrt(t))
}

# Lebesgue density of Z_1
f.z <- function(z, delta, n){
  dnorm(z, mean = sqrt(n) * delta, sd = 1)
  # For t-distribution use
  # dt(z, df = n - 1, ncp = sqrt(n) * delta)
}

# Distribution function of Z_1
F.z <- function(z, delta, n) {
  pnorm(z, mean = sqrt(n) * delta, sd = 1)
  # For t -distribution use
  #pt(z, df = n - 1, ncp = sqrt(n) * delta)
}

# Integration parameters
N <- 5 # Note that N-1 has to be a multiple of 4
wei <- rep(c(32, 12, 32, 14), length.out = N - 1)
wei <- c(7, wei[-length(wei)], 7)
h <- 2 / 45 / (N - 1)


# Introduce class 'design'
design <- function(
  cf,
  ce,
  n1,
  n2,
  c2
) {
  if (!is.numeric(n1)) {
    stop("n1 must be integer")
  }
  return(
    structure(list(
      n1 = n1,
      cf = cf,
      ce = ce,
      n2  = n2,
      c2  = c2
    ), class = c("Design"))
  )
}
