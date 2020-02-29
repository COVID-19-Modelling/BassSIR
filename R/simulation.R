#' Simulate from an estimated Bass-SIR mode
#'
#' @param est a fitted model
#' @param nsim number of simulation
#' @param nforward time-steps forward
#' @param seed random seed
#'
#' @return
#' @export
#'
#' @examples
simulate.estBassSIR <- function(est, nsim = nrow(object$Parameters),
                                nforward = 50, seed = 1167) {
  set.seed(seed)

  n_iter <- nrow(est$Parameters)

  if (nsim == n_iter) {
    selected <- 1:n_iter
  } else if (nsim < n_iter) {
    selected <- sample(1:n_iter, nsim, rep = F)
  } else {
    selected <- sample(1:n_iter, nsim, rep = T)
  }

  f <- switch(est$ModelType,
              BassSIR = system.file("sim/BassSIR.txt", package = "BassSIR"),
              SIR = system.file("sim/SIR.txt", package = "BassSIR"))

  hat <- fitted(est)
  hat <- hat$I_hat[nrow(hat$I_hat), ]

  As <- est$Cases$R + est$Cases$D
  a0 <- As[length(As)]
  dropout <- sum(est$Offsets)

  get_pars <- function(i) {
    pp <- est$Parameters[i, ]

    return(list(
      kappa = pp$kappa,
      beta = pp$beta,
      dropout = dropout,
      m = pp$m,
      I0 = hat[i],
      A0 = a0
    ))
  }

  model <- odin::odin(f)
  cm <- model()

  times <- seq(0, nforward, by = 0.5)

  test <- cm$run(times)
  test <- test[times == round(times), ]

  y0s <- matrix(0, nsim, 2)
  sims <- array(0, c(dim(test), nsim))

  for (i in 1:nsim) {
    key <- selected[i]
    p <- get_pars(key)

    cm$set_user(user = p)
    ys <- cm$run(times)
    ys <- ys[times == round(times), ]
    sims[, , i] <- ys
    y0s[i, ] <- c(p$I0, p$A0)
  }

  colnames(y0s) <- c("I0", "A0")
  rownames(y0s) <- selected

  dimnames(sims)[[2]] <- colnames(test)
  dimnames(sims)[[3]] <- selected

  dates <- est$Cases$Date[est$Cases$len] + 0:nforward

  sim <- list(
    ModelType = est$ModelType,
    Cases = est$Cases,
    Selected = selected,
    Parameters = est$Parameters[selected,],
    Offsets = est$Offsets,
    Y0s = data.frame(y0s),
    Date = dates,
    Simulations = sims
  )

  class(sim) <- "simBassSIR"
  return(sim)
}


#' @rdname simulate.estBassSIR
#' @export
print.simBassSIR <- function(obj) {
  cat("Model type:\t", obj$ModelType, "\n")
  rg <- format.Date(range(obj$Date))
  cat("Time range:\t", rg[1], " - ", rg[2], "\n")
}


#' Summarise epidemiological variables from simulation results
#'
#' @param sim a simulation result
#' @param i an indicator specifying time points, 1 as the default
#'
#' @return
#' @export
#'
#' @examples
summary.simBassSIR <- function(sim, i = 1) {
  loc <- sim$Cases$ID

  indices <- list()
  vs <- sim$Parameters
  indices <- data.table::data.table(Location = loc,
                                    variable = colnames(sim$Parameters),
                                    mean = apply(vs, 2, mean),
                                    lower = apply(vs, 2, quantile, p = 0.025),
                                    upper = apply(vs, 2, quantile, p = 0.975)
  )
  tab <- sim$Simulations

  epis <- matrix(0, nrow(vs), 5)
  colnames(epis) <- c("R0", "R(t)", "PrEx", "PeakSize", "PeakTime")

  epis[, 1] <- tab[i, "R0", ]
  epis[, 2] <- tab[i, "Re", ]
  epis[, 3] <- tab[i, "PrEx", ]
  epis[, 4] <- apply(tab[, "I", ], 2, max)
  epis[, 5] <- apply(tab[, "I", ], 2, which.max)


  indices <- rbind(indices,
                   data.table::data.table(Location = loc,
                                          variable = colnames(epis),
                                          mean = apply(epis, 2, mean),
                                          lower = apply(epis, 2, quantile, p = 0.025),
                                          upper = apply(epis, 2, quantile, p = 0.975)
                   )
  )
  res <- list(
    Location = loc,
    Indices = indices
  )
  class(res) <- "summarySimBassSIR"
  return(res)
}


#' @rdname summary.simBassSIR
#' @export
print.summarySimBassSIR <- function(obj) {
  cat("Location: ", obj$Location, "\n")
  print(obj$Indices[, -1])
}

