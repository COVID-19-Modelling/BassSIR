% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate.estBassSIR}
\alias{simulate.estBassSIR}
\alias{print.simBassSIR}
\alias{aggregate_simulations}
\title{Simulate from an estimated Bass-SIR mode}
\usage{
\method{simulate}{estBassSIR}(est, nsim = nrow(object$Parameters),
  nforward = 50, seed = 1167)

\method{print}{simBassSIR}(obj)

aggregate_simulations(sims)
}
\arguments{
\item{est}{a fitted model}

\item{nsim}{number of simulation}

\item{nforward}{time-steps forward}

\item{seed}{random seed}
}
\description{
Simulate from an estimated Bass-SIR mode
}
