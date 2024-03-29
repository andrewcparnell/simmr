#' Print simmr output object
#'
#' @param x An object of class \code{simmr_output}
#' @param ... Other arguments (not supported)
#'
#' @return Returns a neat summary of the object
#'
#' @seealso \code{\link{simmr_mcmc}} and \code{\link{simmr_ffvb}} for creating
#' \code{simmr_output} objects
#' @export
print.simmr_output <-
  function(x, ...) {
    if (inherits(x, "simmr_output") == TRUE) {
      if (inherits(x, "simmr_mcmc_object") == TRUE) {
        print(x$input)
        message("The input data has been run via simmr_mcmc and has produced ")
        message(nrow(x$output[[1]]$BUGSoutput$sims.matrix), " iterations over ", x$output[[1]]$BUGSoutput$n.chains, " MCMC chains.")
        message("\n\n")
      } else if (inherits(x, "simmr_ffvb_object") == TRUE) {
        print(x$input)
        message("The input data has been run via simmr_ffvb and has produced ")
        message(nrow(x$output[[1]]$BUGSoutput$sims.list$p), " samples.")
      }
    }
  }
