lse <- function(xs, na.rm = F) {
  m <- max(xs)
  return(log(sum(exp(xs - m), na.rm = na.rm)) + m)
}


#' Calculate bootstrap credible interval
#'
#' @param x
#' @param digits
#'
#' @return
#' @export
#'
#' @examples
stats_fn <- function(x, digits = getOption("qwraps2_frmt_digits", 2)) {
  summ <- c(
    mean(x, na.rm = T),
    quantile(x, c(0.025, 0.975), na.rm = T)
  )

  qwraps2::frmtci(summ, digits = digits)
}
