#' Fit a model of the BassSIR family to epidemic data
#'
#' @param d a BassSIR data
#' @param r_rec recovery rate
#' @param r_death death rate
#' @param type type of model
#' @param n_iter number of interation to keep
#' @param ... arguments for jags(...)
#'
#' @return
#' @export
#'
#' @examples
fit <- function(d, r_rec, r_death, type = c("BassSIR", "SIR", "Growth"), hyper, n_iter = 1E4, ...) {
  type <- match.arg(type)

  if (missing(hyper)) hyper <- list()

  # Validate hyperparameters
  hyper$m_mul <- max(hyper$m_mul, 1)
  if (hyper$m_mul <= 1) {
    hyper$m_mul <- 10
  }

  if (type %in% c("BassSIR", "SIR")) {
    hyper$beta_a <- max(hyper$beta_a, 1)
    hyper$beta_b <- max(hyper$beta_b, 1)
  }
  if (type %in% c("BassSIR", "Growth")) {
    hyper$kappa_a <- max(hyper$kappa_a, 1)
    hyper$kappa_b <- max(hyper$kappa_b, 1)
  }


  dat <- (function(nc) {
    ni <- nc$I
    nim2 <- rev(rev(ni)[-(1:2)])
    nim0 <- ni[-(1:2)]
    nr <- (nc$R + nc$D)[-c(1, length(ni))]

    list(
      n_t = length(ni) - 2,
      s = nim0 - nim2,
      x1 = nim0 + nim2,
      x2 = nim0 * nim2,
      A = nr,
      nu = r_rec + r_death,
      mx = max(nc$I + nc$R + nc$D)
    )
  })(d)

  dat <- c(hyper, dat)

  if (type == "Growth") {
    dat$x2 <- NULL
  }

  model.file <- switch(type,
                       BassSIR = system.file("fit/BassSIR.txt", package = "BassSIR"),
                       SIR = system.file("fit/SIR.txt", package = "BassSIR"),
                       Growth = system.file("fit/Growth.txt", package = "BassSIR")
                      )

  inits <- switch(type,
                  BassSIR = function(){ list(tau = 0.01, kappa = 0.01, beta = 0.1, m = dat$mx) },
                  SIR = function(){ list(tau = 0.01, beta = 0.1, m = dat$mx) },
                  Growth = function(){ list(tau = 0.01, kappa = 0.01, m = dat$mx) }
                  )

  pars_save <- c("mu", "m")
  if (type %in% c("BassSIR", "SIR")) pars_save <- c(pars_save, "beta")
  if (type %in% c("BassSIR", "Growth")) pars_save <- c(pars_save, "kappa")

  f <- R2jags::jags(data = dat,
                         inits = inits,
                         parameters.to.save = c("mu", "kappa", "beta", "m"),
                         n.iter = n_iter,
                         model.file = model.file, ...)

  pars_dis <- with(as.data.frame(f$BUGSoutput$sims.matrix), {
    if (type == "Growth") beta <- 0
    if (type == "SIR") kappa <- 0
    data.frame(kappa = kappa, beta = beta, m = m, deviance = deviance)
  })

  pars_con <- with(as.data.frame(f$BUGSoutput$sims.matrix), {
    if (type == "Growth") beta <- 0
    if (type == "SIR") kappa <- 0
    data.frame(kappa = qexp(kappa), beta = qexp(beta), m = m, deviance = deviance)
  })

  res <- list(
    ModelType = type,
    Hyperpars = hyper,
    Parameters = pars_con,
    ParametersDis = pars_dis,
    Cases = cases,
    mus = f$BUGSoutput$sims.matrix[, paste0("mu[", 1:dat$n_t, "]")],
    Offsets = c(r_rec = r_rec, r_death = r_death),
    DIC = f$BUGSoutput$DIC
  )

  class(res) <- paste0("est", type)
  return(res)
}



#' Summarise a fitted Bass SIR model
#'
#' @param est
#'
#' @return
#' @export
#'
#' @examples
summary.estBassSIR <- function(est) {
  y <- list()

  y$Cases <- est$Cases
  y$DIC <- est$DIC
  y$ModelType <- est$ModelType
  y$Pars = with(est$Parameters, {
    list(
      "Kappa*100" = stats_fn(kappa * 100),
      Beta = stats_fn(beta),
      "Effective N (thousand)" = stats_fn(m / 1E3),
      Deviance = stats_fn(deviance))
  })

  class(y) <- "summaryEstBassSIR"
  return(y)
}


#' @rdname summary.estBassSIR
#' @export
print.summaryEstBassSIR <- function(est) {
  print(est$Cases)
  cat("Model: ", est$ModelType, "\n")
  for (name in names(est$Pars)) {
    cat(name, ": \t", est$Pars[[name]], "\n")
  }
}


#' Compare empirical models
#'
#' @param ... models with argument names indicating the names of models
#'
#' @return
#' @export
#'
#' @examples
compare_models <- function(...) {
  m.list <- list(...)

  mods <- names(m.list)
  if (length(mods) < length(m.list)) {
    mods <- paste0("Model_", 1:length(m.list))
  } else {
    mods <- ifelse(mods == "", paste0("Model_", 1:length(m.list)), mods)
  }


  devs <- sapply(m.list, function(x) {
    lse(x$Parameters$deviance / -2) - log(nrow(x$Parameters))
  })

  dics <- sapply(m.list, function(x) x$DIC)

  bfs <- exp(devs - devs[1])

  return(data.frame(
    Location = m.list[[1]]$Cases$ID,
    Models = mods,
    DIC = dics,
    BayesFactor = bfs
  ))
}


