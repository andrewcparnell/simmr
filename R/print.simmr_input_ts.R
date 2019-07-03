#' Print simmr input time series object
#'
#' @param x An object of class \code{simmr_input_ts}
#' @param ... Other arguments (not supported)
#'
#' @return A neat presentation of your simmr object. 
#' @export
print.simmr_input_ts <-
function(x,...) {
  print.simmr_input(x)
  time = x$time_var
  cat('There is a time series variable ranging from value', min(time),'to',max(time),'\n')
  cat('Calling simmr_mcmc wil run a time series simmr model.\n')
}