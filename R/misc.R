lse <- function(xs, na.rm = F) {
  m <- max(xs)
  return(log(sum(exp(xs - m), na.rm = na.rm)) + m)
}


#' Calculate bootstrap credible interval
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
stats_fn <- function(x) {
  summ <- c(
    mean(x, na.rm = T),
    quantile(x, c(0.025, 0.975), na.rm = T)
  )

  qwraps2::frmtci(summ)
}


calc_bayes_factor <- function(...) {
  m.list <- list(...)
  devs <- sapply(m.list, function(x) {
    lse(x$Parameters$deviance / -2) - log(nrow(x$Parameters))
  })

  dics <- sapply(m.list, function(x) x$DIC)
  names(dics) <- sapply(m.list, function(x) paste0("DIC.", x$Model))

  bfs <- exp(devs - devs[1])
  names(bfs) <- sapply(m.list, function(x) paste0("BF.", x$Model))

  return(c(
    list(Location = m.list[[1]]$Summary$Location),
    as.list(bfs),
    as.list(dics)
  ))
}
