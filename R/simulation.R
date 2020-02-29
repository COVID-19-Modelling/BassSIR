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
  cat("Time ragne:\t", rg[1], " - ", rg[2], "\n")
}
