#' Fit a model of the BassSIR family to epidemic data
#'
#' @param d a BassSIR data
#' @param r_rec recovery rate
#' @param r_death death rate
#' @param type type of model
#' @param n_iter number of interation
#' @param n_collect number of interation to collect
#' @param ... arguments for jags(...)
#'
#' @return
#' @export
#'
#' @examples
fit <- function(d, r_rec, r_death, type = c("BassSIR", "SIR", "Growth"), hyper, n_iter = 1E4, n_collect = 1000, ...) {
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
      di = nim0 - nim2,
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
                    parameters.to.save = pars_save,
                    n.iter = n_iter,
                    n.burnin = n_iter - n_collect,
                    model.file = model.file, ...)

  pars_dis <- with(as.data.frame(f$BUGSoutput$sims.matrix), {
    if (type == "Growth") beta <- 0
    if (type == "SIR") kappa <- 0
    data.frame(kappa = kappa, beta = beta, m = m, r_death = r_death, r_rec = r_rec, deviance = deviance)
  })

  pars_con <- with(as.data.frame(f$BUGSoutput$sims.matrix), {
    if (type == "Growth") beta <- 0
    if (type == "SIR") kappa <- 0
    data.frame(kappa = qexp(kappa), beta = qexp(beta), m = m, r_death = r_death, r_rec = r_rec, deviance = deviance)
  })

  res <- list(
    ModelType = type,
    Hyperpars = hyper,
    Parameters = pars_con,
    ParametersDis = pars_dis,
    Cases = d,
    mus = f$BUGSoutput$sims.matrix[, paste0("mu[", 1:dat$n_t, "]")],
    Offsets = c(r_rec = r_rec, r_death = r_death),
    BUGS = f$BUGSoutput,
    DIC = f$BUGSoutput$DIC
  )

  class(res) <- "estBassSIR"
  return(res)
}


#' Fit a full model of the BassSIR family to epidemic data
#'
#' @param d a BassSIR data
#' @param type type of model
#' @param r_iso rate to isolation
#' @param n_iter number of interation
#' @param n_collect number of interation to collect
#' @param ... arguments for jags(...)
#'
#' @return
#' @export
#'
#' @examples
fit_full <- function(d, r_iso = 1/6.1, type = c("BassSIR", "SIR", "Growth"), hyper, n_iter = 1E4, n_collect = 1000, ...) {
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
    nh <- nc$I
    nhm2 <- rev(rev(nh)[-(1:2)])
    nhm0 <- nh[-(1:2)]

    nd <- nc$D
    ndm2 <- rev(rev(nd)[-(1:2)])
    ndm0 <- nd[-(1:2)]

    nr <- nc$R
    nrm2 <- rev(rev(nr)[-(1:2)])
    nrm0 <- nr[-(1:2)]

    list(
      n_t = length(nh) - 2,
      dh = nhm0 - nhm2,
      dd = ndm0 - ndm2,
      dr = nrm0 - nrm2,
      x1 = nhm0 + nhm2,
      # x2 = nhm0 * nhm2,
      A = (nc$R + nc$D + nc$I)[-c(1, length(nh))],
      r_iso = r_iso,
      mx = max(nc$I + nc$R + nc$D)
    )
  })(d)


  dat <- c(hyper, dat)

  if (type == "Growth") {
    dat$x2 <- NULL
  }

  model.file <- switch(type,
                       BassSIR = system.file("fit/FullBassSIR.txt", package = "BassSIR"),
                       SIR = system.file("fit/FullSIR.txt", package = "BassSIR"),
                       Growth = system.file("fit/FullGrowth.txt", package = "BassSIR")
  )

  inits <- switch(type,
                  BassSIR = function(){ list(tau_h = 0.01, kappa = 0.01, beta = 0.1, m = dat$mx) },
                  SIR = function(){ list(tau_h = 0.01, beta = 0.1, m = dat$mx) },
                  Growth = function(){ list(tau_h = 0.01, kappa = 0.01, m = dat$mx) }
  )

  pars_save <- c("I", "mu", "m", "r_iso", "r_death", "r_rec")

  if (type %in% c("BassSIR", "SIR")) pars_save <- c(pars_save, "beta")
  if (type %in% c("BassSIR", "Growth")) pars_save <- c(pars_save, "kappa")

  f <- R2jags::jags(data = dat,
                    inits = inits,
                    parameters.to.save = pars_save,
                    n.iter = n_iter,
                    n.burnin = n_iter - n_collect,
                    model.file = model.file, ...)

  pars_dis <- with(as.data.frame(f$BUGSoutput$sims.matrix), {
    if (type == "Growth") beta <- 0
    if (type == "SIR") kappa <- 0
    data.frame(kappa = kappa, beta = beta, m = m,
               r_death = r_death, r_rec = r_rec, r_iso = r_iso, deviance = deviance)
  })

  pars_con <- with(as.data.frame(f$BUGSoutput$sims.matrix), {
    if (type == "Growth") beta <- 0
    if (type == "SIR") kappa <- 0
    data.frame(kappa = qexp(kappa), beta = qexp(beta), m = m,
               r_death = r_death, r_rec = r_rec, r_iso = r_iso, deviance = deviance)
  })

  res <- list(
    ModelType = type,
    Hyperpars = hyper,
    Parameters = pars_con,
    ParametersDis = pars_dis,
    Cases = d,
    mus = f$BUGSoutput$sims.matrix[, paste0("mu[", 1:dat$n_t, "]")],
    is = f$BUGSoutput$sims.matrix[, paste0("I[", 1:dat$n_t + 1, "]")],
    BUGS = f$BUGSoutput,
    DIC = f$BUGSoutput$DIC
  )

  class(res) <- "estFullBassSIR"
  return(res)
}



#' Summarise a fitted Bass SIR model
#'
#' @param est a model with fitted parameters
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
    r0 <- beta / (r_death + r_rec)
    rt <- r0 * (1 - c(est$Cases$I + est$Cases$R + est$Cases$D)[est$Cases$len] / m)

    list(
      "Kappa*100" = stats_fn(kappa * 100),
      Beta = stats_fn(beta),
      "Effective N (thousand)" = stats_fn(m / 1E3),
      R0 = stats_fn(r0),
      "R(t)" = stats_fn(rt),
      DeathRate = stats_fn(r_death, 4),
      RecoveryRate = stats_fn(r_rec, 4),
      Deviance = stats_fn(deviance))
  })

  class(y) <- "summaryEstBassSIR"
  return(y)
}


#' @rdname summary.estBassSIR
#' @export
summary.estFullBassSIR <- function(est) {
  y <- list()

  y$Cases <- est$Cases
  y$DIC <- est$DIC
  y$ModelType <- est$ModelType

  y$Pars = with(est$Parameters, {
    r0 <- beta / r_iso
    rt <- r0 * (1 - (c(est$Cases$I + est$Cases$R + est$Cases$D)[est$Cases$len] + est$is[, ncol(est$is)]) / m)

    list(
      "Kappa*100" = stats_fn(kappa * 100),
      Beta = stats_fn(beta),
      "Effective N (thousand)" = stats_fn(m / 1E3),
      R0 = stats_fn(r0),
      "R(t)" = stats_fn(rt),
      DeathRate = stats_fn(r_death, 4),
      RecoveryRate = stats_fn(r_rec, 4),
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
    mods <- ifelse(mods == "" | is.null(mods), paste0("Model_", 1:length(m.list)), mods)
  }

  deviance = sapply(m.list, function(x) {
    stats_fn(x$Parameters$deviance)
  })

  devs <- sapply(m.list, function(x) {
    lse(x$Parameters$deviance / -2) - log(nrow(x$Parameters))
  })

  dics <- sapply(m.list, function(x) x$DIC)

  bfs <- exp(devs - devs[1])

  return(data.table::data.table(
    Location = m.list[[1]]$Cases$ID,
    Model = mods,
    Deviance = deviance,
    DIC = dics,
    BayesFactor = bfs
  ))
}


