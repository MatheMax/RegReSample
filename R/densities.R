# Prior density
pi_0 <- function(delta, delta.alt){
  dnorm(delta, mean = delta.alt, sd = 0.1)
}

# Posterior density
pi_1 <- function(delta, delta.alt){
  dnorm(delta, mean = delta.alt, sd = 0.1)
}

# Lebesgue density of Z_1
f.z <- function(z, delta, n){
  dnorm(z, mean = sqrt(n) * delta, sd = 1)
}

# Distribution function of Z_1
F.z <- function(z, delta, n) {
  pnorm(z, mean = sqrt(n) * delta, sd = 1)
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
