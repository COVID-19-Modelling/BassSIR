#' Fitted values of BassSIR models
#'
#' @param est
#' @param method
#'
#' @return
#' @export
#'
#' @examples
fitted.estBassSIR <- function(est, method = c("match", "forward")) {
  mus <- t(est$mus)

  method <- match.arg(method)

  nc <- length(est$Cases$I)

  if (method == "forward") {
    hat <- est$Cases$I[1:(nc - 2)] + mus
  } else {
    hat <- (est$Cases$I[1:(nc - 2)] + est$Cases$I[3:nc] + mus) / 2
  }

  r0 <- est$Parameters$beta / sum(est$Offsets)
  a <- (est$Cases$R + est$Cases$D)[-c(1:2)]
  rt <- r0 * (1 - (hat + a) / est$Parameters$m)

  y <- list(
    I_data = est$Cases$I[1:(nc - 2)],
    I_hat = hat,
    Rt_hat = rt,
    method = method
  )

  class(y) <- "fittedBassSIR"
  return(y)
}


#' @rdname fitted.EstBassSIR
#' @export
summary.fittedBassSIR <- function(f) {

}
