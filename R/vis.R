#' Combine information for visualise prevalence
#'
#' @param est model estimates
#' @param sim model simulations
#' @param method method for finding fitted values
#'
#' @return
#' @export
#'
#' @examples
to_vis_trajectory <- function(est, sim, traj = c("I", "Rt"), method = "match") {
  f <- fitted(est, method = method)
  traj <- match.arg(traj)

  if (traj == "I") {
    traj_back <- data.table::data.table(Time = sim$Date[1] - nrow(f$I_hat):1 + 1,
                                 mean = apply(f$I_hat, 1, mean),
                                 lower = apply(f$I_hat, 1, quantile, p = 0.025),
                                 upper = apply(f$I_hat, 1, quantile, p = 0.975))

    traj_fore <- data.table::data.table(Time = sim$Date,
                                        mean = apply(sim$Simulations[, "I",], 1, mean),
                                        lower = apply(sim$Simulations[, "I",], 1, quantile, p = 0.025),
                                        upper = apply(sim$Simulations[, "I",], 1, quantile, p = 0.975))

    dat <- data.table::data.table(Time = sim$Date[1] - length(f$I_data):1 + 1, Data = f$I_data)

  } else {
    traj_back <- data.table::data.table(Time = sim$Date[1] - nrow(f$I_hat):1 + 1,
                                        mean = apply(f$Rt_hat, 1, mean),
                                        lower = apply(f$Rt_hat, 1, quantile, p = 0.025),
                                        upper = apply(f$Rt_hat, 1, quantile, p = 0.975))

    traj_fore <- data.table::data.table(Time = sim$Date,
                                        mean = apply(sim$Simulations[, "Re",], 1, mean),
                                        lower = apply(sim$Simulations[, "Re",], 1, quantile, p = 0.025),
                                        upper = apply(sim$Simulations[, "Re",], 1, quantile, p = 0.975))

    dat <- NULL
  }

  list(
    Data = dat,
    Fitted = traj_back,
    Simulation = traj_fore
  )
}
