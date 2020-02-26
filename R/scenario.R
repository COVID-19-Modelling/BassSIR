run_scenario <- function(sim, fn) {

  f <- switch(est$ModelType,
              BassSIR = system.file("sim/IntvBassSIR.txt", package = "BassSIR"),
              SIR = system.file("sim/IntvSIR.txt", package = "BassSIR"))

  y0s <- sim$Y0s
  pars <- sim$Parameters


  tt <- 0:length(sim$Date)
  dropout <- est$Offsets)

  get_pars <- function(i) {
    pp <- pars[i, ]
    y0 <- y0s[i, ]

    return(list(
      kappa = pp$kappa,
      beta = pp$beta,
      dropout = dropout,
      m = pp$m,
      I0 = hat[i],
      A0 = a0
    ))
  }

  times <- seq(0, length(sim$Date), 0.5)
  model <- odin::odin(f)
  cm <- model()



}


compare_scenarios <- function(sim, ...) {
  scs <- lsit(...)
}
