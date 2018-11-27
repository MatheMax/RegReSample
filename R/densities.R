'
f.z <- function(z, delta, n){
  dnorm(z, mean = sqrt(n) * delta, sd = 1)  # Lebesgue density of Z_1
}
F.z <- function(z, delta, n) {
  pnorm(z, mean = sqrt(n) * delta, sd = 1)  # Distribution function of Z_1
}
pi.0 <- function(delta, delta.alt, tau){
  dnorm(delta, mean = delta.alt, sd = tau)
}
pi.1 <- function(delta, delta.alt, tau, z1, n1){
  dnorm(delta, mean = delta.alt, sd = tau)
}
'


# Integration parameters
N <- 5 # Note that N-1 has to be a multiple of 4
wei <- rep(c(32, 12, 32, 14), length.out = N - 1)
wei <- c(7, wei[-length(wei)], 7)
h <- 2 / 45 / (N - 1)

N.2 <- 41
wei.2 <- rep(c(32, 12, 32, 14), length.out = N.2 - 1)
wei.2 <- c(7, wei.2[-length(wei.2)], 7)
h.2 <- 2 / 45 / (N.2 - 1)

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
