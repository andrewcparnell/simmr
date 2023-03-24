#' Summarises the output created with \code{\link{simmr_mcmc_tdf}}
#'
#' Produces textual summaries and convergence diagnostics for an object created
#' with \code{\link{simmr_mcmc_tdf}}. The different options are: 'diagnostics'
#' which produces Brooks-Gelman-Rubin diagnostics to assess MCMC convergence,
#' 'quantiles' which produces credible intervals for the parameters,
#' 'statistics' which produces means and standard deviations, and
#' 'correlations' which produces correlations between the parameters.
#'
#' The quantile output allows easy calculation of 95 per cent credible
#' intervals of the posterior dietary proportions. The Gelman diagnostic values should be close to 1 to ensure satisfactory convergence.
#'
#' Multiple groups are not currently supported for estimating TDFs.
#'
#' @param object An object of class \code{simmr_output_df} produced by the
#' function \code{\link{simmr_mcmc_tdf}}
#' @param type The type of output required. At least none of 'diagnostics',
#' 'quantiles', 'statistics', 'correlations'
#' @param ...  Not used
#' @return An list containing the following components: \item{gelman }{The
#' convergence diagnostics} \item{quantiles }{The quantiles of each parameter
#' from the posterior distribution} \item{statistics }{The means and standard
#' deviations of each parameter} \item{correlations }{The posterior
#' correlations between the parameters} Note that this object is reported
#' silently so will be discarded unless the function is called with an object
#' as in the example below.
#'
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#'
#' @seealso See \code{\link{simmr_mcmc_tdf}} for creating objects suitable for this
#' function, and many more examples. See also \code{\link{simmr_load}} for
#' creating simmr objects, \code{\link{plot.simmr_input}} for creating isospace
#' plots, \code{\link{plot.simmr_output}} for plotting output.
#'
#' @importFrom stats sd cor
#'
#' @export
summary.simmr_output_tdf <-
  function(object, type = c("diagnostics", "quantiles", "statistics", "correlations"), ...) {
    # Get the specified type
    type <- match.arg(type, several.ok = TRUE)

    out_all <- object$output$BUGSoutput$sims.matrix

    # Get objects
    out_bgr <- object$output$BUGSoutput$summary[, "Rhat"]
    out_quantiles <- t(apply(out_all, 2, "quantile", probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    #  coda:::summary.mcmc.list(object$output)$quantiles
    out_statistics <- t(apply(out_all, 2, function(x) {
      return(c(mean = mean(x), sd = stats::sd(x)))
    }))
    # coda:::summary.mcmc.list(object$output)$statistics[,1:2]
    out_cor <- stats::cor(out_all)

    if ("diagnostics" %in% type) {
      # Print out gelman diagnostics of the output
      cat("Gelman diagnostics - these values should all be close to 1.\n")
      cat("If not, try a longer run of simmr_mcmc_tdf.\n")
      print(round(out_bgr, 2))
    }

    if ("quantiles" %in% type) {
      # Print out quantiles argument
      print(round(out_quantiles, 3))
    }

    if ("statistics" %in% type) {
      # Print out quantiles argument
      print(round(out_statistics, 3))
    }

    if ("correlations" %in% type) {
      # Print out quantiles argument
      print(round(out_cor, 3))
    }

    invisible(list(gelman = out_bgr, quantiles = out_quantiles, statistics = out_statistics, correlations = out_cor))
  }
