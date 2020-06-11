#' Print a simmr output object
#'
#' @param x An object of class \code{simmr_output}
#' @param ... Other arguments (not supported)
#'
#' @return Returns a neat summary of the object
#' 
#' @seealso \code{\link{simmr_mcmc}} for creating \code{simmr_output} objects
#' @export
print.simmr_output <-
function(x,...) {
  print(x$input)
  cat('The input data has been run via simmr_mcmc and has produced ')
  cat(nrow(x$output$BUGSoutput$sims.matrix),'iterations over',x$output$BUGSoutput$n.chains,'MCMC chains.')
  cat('\n\n')
}
