#' Plot important key figures of a two-stage design
#'
#' Returns 6 plots for the first-stage,
#' n_2-function, c_2-function,
#' power-curve, conditional power curve
#' and expected sample size
#'
#' @param d An object of class design
#'
#' @export

plot_design <- function(d, distribution = c("normal", "t"), delta.alt = .3,
                        prior = c("point", "normal"),
                        power = c("predictive", "equal"),
                        delta.mcr = .1,
                        beta = .2,
                        beta.2 = .2,
                        tau = .1) {
  library(dplyr)
  library(ggplot2)

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
    pi.1 <- function(delta, z1){
      dnorm(delta, mean = delta.alt, sd = tau)
    }
  } else if(power == "predictive"){
    pi.1 <- function(delta, z1){
      t <- 1 / (d$n1 + 1 / tau^2)
      dnorm(delta, mean = t * (delta.alt / tau^2 + sqrt(d$n1) * z1), sd = sqrt(t))
    }
  }


  # First stage
  df <- dplyr::data_frame(
    ` ` = c("n[1]", "c[f]", "c[e]"),
    `Value` = c(2*ceiling(d$n1), round(d$cf, 2), round(d$ce, 2))
  )

  table <- gridExtra::tableGrob(df, rows = NULL,
                                theme=gridExtra::ttheme_minimal(base_size=11, parse=TRUE))
  title <- grid::textGrob("First stage", gp = grid::gpar(fontsize=12))
  padding <- grid::unit(0.8, "line")
  table <- gtable::gtable_add_rows(
    table, heights = grid::grobHeight(title) + padding, pos = 0
  )
  table <- gtable::gtable_add_grob(
    table, list(title),
    t = 1, l = 1, r = ncol(table)
  )
  bild1 <- table


  z <- seq(round(d$cf, 2) - 0.2, round(d$ce, 2) + 0.2, 0.1)

  z1 <- seq(min(z), d$cf, 0.1)
  z2 <- seq(d$cf, d$ce, 0.1)
  z3 <- seq(d$ce, max(z), 0.1)


  # Stage two sample size
  v1 <- rep(0, length(z1))
  v2 <- sapply(z2, function(x){2 * (d$n2(x))})
  v3 <- rep(0, length(z3))


  out2 <- suppressWarnings(data.frame(cbind(z1, z2, z3, v1, v2, v3)))
  bild2 <-   ggplot(out2) +
    geom_line(aes(x = z1, y = v1, linetype = "solid"), size = 0.3) +
    geom_line(aes(x = z2, y = v2, linetype = "solid"), size = 0.3) +
    geom_line(aes(x = z3, y = v3, linetype = "solid"), size = 0.3) +
    ggplot2::labs(title = expression(n[2](z[1])),
                  x = expression(z[1]),
                  y = element_blank()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = element_text(size=12),
      legend.position = "none"
    ) +
    scale_x_continuous(breaks = seq(-2, 4, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1000, 50)) +
    expand_limits(y = c(0, max(v2)))


  # Stage two rejection boundary
  v22 <- sapply(z2, function(x){d$c2(x)})

  out3 <- suppressWarnings(data.frame(cbind(z2, v22)))
  bild3 <-   ggplot(out3) +
    geom_line(aes(x = z2, y = v22, linetype = "solid"), size = 0.3) +
    ggplot2::labs(title = expression(c[2](z[1])),
                  x = expression(z[1]),
                  y=element_blank()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = element_text(size=12),
      legend.position = "none"
    ) +
    scale_x_continuous(breaks = seq(-2, 4, 0.5)) +
    scale_y_continuous(breaks = seq(-2, 4, 0.5)) +
    expand_limits(y = c(floor(min(v22)), ceiling(max(v22))))

  # Power
  theta <- seq(0, max(.5, delta.alt+.05), 0.05)
  pow.plot <- function(th) {
    1 - F.z(d$cf, th, d$n1) -
      integrate(Vectorize(function(z){
        F.z(d$c2(z), th, d$n2(z)) * f.z(z, th, d$n1)
      }),
      d$cf,
      d$ce)$value
  }
  v41 <- sapply(theta, function(x){pow.plot(x)})
  out4 <- data.frame(cbind(theta, v41))
  bild4 <-   ggplot2::ggplot(out4, ggplot2::aes(x = theta)) +
    ggplot2::geom_line(ggplot2::aes(y = v41, linetype = "solid"), size=0.3) +
    ggplot2::geom_hline(yintercept = 1-beta, colour = "grey", size=0.3) +
    ggplot2::geom_vline(xintercept = delta.alt, colour = "grey", size=0.3) +
    ggplot2::labs(title = "Global power for point alternative",
                  x = expression(delta),
                  y = element_blank()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = element_text(size=12),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.1))


  # Expected sample size
  plot.expn <- function(th) {
    d$n1 + integrate(Vectorize(function(z){
      d$n2(z) * f.z(z, th, d$n1)
    }),
    d$cf,
    d$ce)$value
  }
  v51 <- sapply(theta, function(x){2 * plot.expn(x)})
  out5 <- data.frame(cbind(theta, v51))
  bild5 <-   ggplot2::ggplot(out5, ggplot2::aes(x = theta)) +
    ggplot2::geom_line(ggplot2::aes(y = v51, linetype = "solid"), size=0.3) +
    ggplot2::geom_vline(xintercept = delta.alt, colour = "grey", size=0.3) +
    ggplot2::labs(title = "Expected sample size",
                  x = expression(delta),
                  y = element_blank()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = element_text(size=12),
      legend.position = "none"
    ) +
    expand_limits(y = c(floor(min(v51)), ceiling(max(v51)))) +
    scale_y_continuous(breaks = seq(0, 1000, 10))


  # Conditional power

  plot.cp <- function(z, delta){
    pnorm(sqrt(d$n2(z)) * delta - d$c2(z))
  }

  v61 <- rep(0, length(z1))
  if(prior == "normal"){
    if(power == "predictive"){
      v62 <- sapply(z2, function(x){
        integrate(Vectorize(function(delta){plot.cp(x, delta) * pi.1(delta, x)}),
                  delta.mcr,
                  Inf)$value
      })
    } else if(power == "equal"){
      v62 <- sapply(z2, function(x){
        integrate(Vectorize(function(delta){plot.cp(x, delta) * pi.0(delta)}),
                  delta.mcr,
                  Inf)$value
    })
    }
  } else{
    v62 <- sapply(z2, function(x){plot.cp(x, delta.alt)})
  }
  v63 <- rep(1, length(z3))

  out6 <- suppressWarnings(data.frame(cbind(z1, z2, z3, v61, v62, v63)))
  bild6 <- ggplot(out6) +
    ggplot2::geom_hline(yintercept = 1-beta.2, colour = "grey", size=0.3) +
    ggplot2::geom_line(aes(x = z1, y = v61, linetype = "solid"), size=0.3) +
    ggplot2::geom_line(aes(x = z2, y = v62, linetype = "solid"), size=0.3) +
    ggplot2::geom_line(aes(x = z3, y = v63, linetype = "solid"), size=0.3) +
    ggplot2::labs(title = "Conditional power",
                  x = expression(z[1]),
                  y = element_blank()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = element_text(size=12),
      legend.position = "none"
    ) +
    scale_x_continuous(breaks = seq(-2, 4, 0.5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1))


  return(gridExtra::grid.arrange(bild1, bild2, bild3, bild4, bild5, bild6, ncol=3))

}
