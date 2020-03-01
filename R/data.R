#' Confirmed cases, the recovered, and peopole who died from COVID-19
#'
#' A dataset containing the Confirmed cases, the recovered, and peopole who died from COVID-19
#' by chinese province
#'
#'
#' @format A list of data frames Province -> data.frame:
#' \describe{
#'   \item{Date}{Date of the entry}
#'   \item{Confirmed}{Confirmed cases}
#'   \item{Dead}{People who died from COVID-19}
#'   \item{Cured}{People who recovered from COVID-19}
#'   ...
#' }
#' @source \url{}
"n_covid19"



#' Prepare data for BassSIR modelling
#'
#' @param d source data.frame
#' @param location identity of the data
#' @param status the source data are "Cumulative" numbers or "Current"
#' @param i_date index indicating date
#' @param i_confirmed index indicating confirmed cases
#' @param i_recovered index indicating recovered people
#' @param i_dead index indicating the number of dead people
#'
#' @return
#' @export
#'
#' @examples
#' as_bass_data(n_covid19$Hubei)
#'
as_bass_data <- function(d, id = "Place X", status = c("Current", "Cumulative"), i_date = "Date",
                         i_confirmed = "Confirmed", i_recovered = "Cured", i_dead = "Dead") {

  if (!is.data.frame(d)) {
    d <- as.data.frame(d)
  }
  d <- d[order(d[i_date]), ]


  dates <- range(d[, i_date])

  ind_date <- d[, i_date] - dates[1]
  dates <- dates[1] + 0:diff(dates)

  xout <- 1:length(dates) - 1

  status <- match.arg(status)
  cases <- list(
    ID = id,
    Date = dates,
    len = length(dates),
    I = ts(approx(ind_date, d[, i_confirmed], xout, method = "const")$y),
    R = ts(approx(ind_date, d[, i_recovered], xout, method = "const")$y),
    D = ts(approx(ind_date, d[, i_dead], xout, method = "const")$y)
  )

  if (status == "Cumulative") {
    cases$I <- cases$I - cases$R - cases$D
  }

  class(cases) <- "BassData"
  return(cases)
}


#' @rdname as_bass_data
#' @export
print.BassData <- function(cases) {
  cat("Time-series of ", cases$ID, "\n")
  cat("From: ", as.character(cases$Date[1]),
      ", to: ", as.character(cases$Date[cases$len]), "\n")
}
