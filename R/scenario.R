valid_pars <- function(pars) {
  nt <- length(pars$tt)
  all(sapply(pars, length) == nt)
}


#' Run scenario analysis
#'
#' @param sim baseline simulations
#' @param fn a function for modifying parameters
#'
#' @return
#' @export
#'
#' @examples
#' library(R2jags)
#' cases <- as_bass_data(n_covid19$Hubei)
#' est <- BassSIR::fit(cases, r_rec = 1/22.2, r_death = 1/22.3, type = "BassSIR")
#'
#' sim <- simulate(est, nsim = 500)
#'
#' zero_kappa <- function(pars) {
#'   pars$kappa <- rep(0, length(pars$kappa))
#'   return(pars)
#' }
#'
#' lockdown <- run_scenario(sim, zero_kappa)
#'
#' cp <- compare_scenarios(sim, Lockdown = lockdown)
run_scenario <- function(sim, fn) {

  f <- switch(sim$ModelType,
              BassSIR = system.file("sim/IntvBassSIR.txt", package = "BassSIR"),
              SIR = system.file("sim/IntvSIR.txt", package = "BassSIR"))

  y0s <- sim$Y0s

  get_pars <- function(i) {
    p <- sim$Parameters[i, ]
    nt <- length(sim$Date)

    list(
      tt = 1:nt - 1,
      m = rep(p$m, nt),
      beta = rep(p$beta, nt),
      kappa = rep(p$kappa, nt),
      r_rec = rep(as.numeric(sim$Offsets[1]), nt),
      r_death = rep(as.numeric(sim$Offsets[2]), nt)
    )
  }

  pars <- get_pars(1)
  times <- seq(min(pars$tt), max(pars$tt), 0.5)

  pars <- fn(pars)
  stopifnot(valid_pars(pars))

  model <- odin::odin(f)
  cm <- model(user = pars)
  test <- cm$run(times)
  test <- test[times == round(times), ]

  sims <- array(0, c(dim(test), nrow(sim$Parameters)))

  for (i in 1:nrow(sim$Parameters)) {
    pars <- fn(get_pars(i))
    pars$I0 <- sim$Y0s$I0[i]
    pars$A0 <- sim$Y0s$A0[i]

    cm$set_user(user = pars)
    ys <- cm$run(times)
    ys <- ys[times == round(times), ]
    sims[, , i] <- ys
  }

  dimnames(sims)[[2]] <- colnames(test)
  dimnames(sims)[[3]] <- dimnames(sim$Simulations)[[3]]

  res <- list(
    Simulations = sims,
    Date = sim$Date
  )

  class(res) <- "simScBassSIR"
  return(res)
}


#' @rdname run_scenario
#' @export
compare_scenarios <- function(sim, ..., fn_change = c("Y0", "Y1", "Yt")) {
  scs <- list(...)


  scs <- list(lockdown = lockdown)

  vs <- intersect(dimnames(sim$Simulations)[[2]], dimnames(scs[[1]]$Simulations)[[2]])
  vs <- vs[vs != "t"]

  traj <- list()

  for (v in vs) {
    sims <- sim$Simulations[, v, ]
    temp <- data.frame(Time = sim$Date, Scenario = "Baseline",
                       mean = apply(sims, 1, mean),
                       lower = apply(sims, 1, quantile, p = 0.025),
                       upper = apply(sims, 1, quantile, p = 0.975), row.names = NULL)
    for (sc in names(scs)) {
      sims <- scs[[sc]]$Simulations[, v, ]
      temp <- rbind(temp,
                    data.frame(Time = sim$Date, Scenario = sc,
                               mean = apply(sims, 1, mean),
                               lower = apply(sims, 1, quantile, p = 0.025),
                               upper = apply(sims, 1, quantile, p = 0.975), row.names = NULL)
      )
    }
    traj[[v]] <- temp
  }


  change <- list()
  fn_change <- match.arg(fn_change)

  for (v in vs) {
    baseline <- sim$Simulations[, v, ]

    temp <- data.frame(Time = sim$Date, Scenario = "Baseline",
                       mean = 0, lower = 0, upper = 0, row.names = NULL)
    for (sc in names(scs)) {
      sims <- scs[[sc]]$Simulations[, v, ]
      sims <- switch (fn_change,
        Y0 = sims / baseline - 1,
        Y1 = 1 - baseline / sims,
        Yt = log(sims) - log(sims)
      ) * 100
      temp <- rbind(temp,
                    data.frame(Time = sim$Date, Scenario = sc,
                               mean = apply(sims, 1, mean),
                               lower = apply(sims, 1, quantile, p = 0.025),
                               upper = apply(sims, 1, quantile, p = 0.975), row.names = NULL)
      )
    }
    change[[v]] <- temp
  }

  res <- list(
    Trajectories = traj,
    Changes = change,
    Scenarios = names(scs)
  )
  class(res) <- "compare_scenarios"
  return(res)
}
