#' Summarises the output created with \code{\link{simmr_mcmc}}
#'
#' Produces textual summaries and convergence diagnostics for an object created
#' with \code{\link{simmr_mcmc}}. The different options are: 'diagnostics'
#' which produces Brooks-Gelman-Rubin diagnostics to assess MCMC convergence,
#' 'quantiles' which produces credible intervals for the parameters,
#' 'statistics' which produces means and standard deviations, and
#' 'correlations' which produces correlations between the parameters.
#'
#' The quantile output allows easy calculation of 95 per cent credible
#' intervals of the posterior dietary proportions. The correlations, along with
#' the matrix plot in \code{\link{plot.simmr_output}} allow the user to judge
#' which sources are non-identifiable. The Gelman diagnostic values should be
#' close to 1 to ensure satisfactory convergence.
#'
#' When multiple groups are included, the output automatically includes the
#' results for all groups.
#'
#' @param object An object of class \code{simmr_output} produced by the
#' function \code{\link{simmr_mcmc}}
#' @param type The type of output required. At least none of 'diagnostics',
#' 'quantiles', 'statistics', or 'correlations'.
#' @param group Which group or groups the output is required for.
#' @param ...  Not used
#' @return A list containing the following components: \item{gelman }{The
#' convergence diagnostics} \item{quantiles }{The quantiles of each parameter
#' from the posterior distribution} \item{statistics }{The means and standard
#' deviations of each parameter} \item{correlations }{The posterior
#' correlations between the parameters} Note that this object is reported
#' silently so will be discarded unless the function is called with an object
#' as in the example below.
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} for creating objects suitable for this
#' function, and many more examples. See also \code{\link{simmr_load}} for
#' creating simmr objects, \code{\link{plot.simmr_input}} for creating isospace
#' plots, \code{\link{plot.simmr_output}} for plotting output.
#'
#' @importFrom stats sd cor
#'
#' @examples
#' \dontrun{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
#' # The data
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' # Plot
#' plot(simmr_1)
#'
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Summarise
#' summary(simmr_1_out) # This outputs all the summaries
#' summary(simmr_1_out, type = "diagnostics") # Just the diagnostics
#' # Store the output in an
#' ans <- summary(simmr_1_out,
#'   type = c("quantiles", "statistics")
#' )
#' }
#' @export
summary.simmr_output <-
  function(object, type = c("diagnostics", "quantiles", "statistics", "correlations"), group = 1, ...) {

    # Get the specified type
    type <- match.arg(type, several.ok = TRUE)

    # Set up containers
    out_bgr <- out_quantiles <- out_statistics <- out_cor <- vector("list", length = length(group))
    group_names <- levels(object$input$group)
    names(out_bgr) <- paste0("group_", group)
    names(out_quantiles) <- paste0("group_", group)
    names(out_statistics) <- paste0("group_", group)
    names(out_cor) <- paste0("group_", group)

    # Loop through groups
    for (i in 1:length(group)) {
      cat(paste("\nSummary for", group_names[group[i]], "\n"))
      out_all <- object$output[[group[i]]]$BUGSoutput$sims.matrix

      # Get objects
      out_bgr[[i]] <- object$output[[i]]$BUGSoutput$summary[, "Rhat"]
      out_quantiles[[i]] <- t(apply(out_all, 2, "quantile", probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
      #  coda:::summary.mcmc.list(object$output)$quantiles
      out_statistics[[i]] <- t(apply(out_all, 2, function(x) {
        return(c(mean = mean(x), sd = stats::sd(x)))
      }))
      # coda:::summary.mcmc.list(object$output)$statistics[,1:2]
      out_cor[[i]] <- stats::cor(out_all)

      if ("diagnostics" %in% type) {
        # Print out gelman diagnostics of the output
        cat("Gelman diagnostics - these values should all be close to 1.\n")
        cat("If not, try a longer run of simmr_mcmc.\n")
        print(round(out_bgr[[i]], 2))
      }

      if ("quantiles" %in% type) {
        # Print out quantiles argument
        print(round(out_quantiles[[i]], 3))
      }

      if ("statistics" %in% type) {
        # Print out quantiles argument
        print(round(out_statistics[[i]], 3))
      }

      if ("correlations" %in% type) {
        # Print out quantiles argument
        print(round(out_cor[[i]], 3))
      }
    }

    if (object$input$n_groups == 1) {
      invisible(list(gelman = out_bgr[[1]], quantiles = out_quantiles[[1]], statistics = out_statistics[[1]], correlations = out_cor[[1]]))
    } else {
      invisible(list(gelman = out_bgr, quantiles = out_quantiles, statistics = out_statistics, correlations = out_cor))
    }
  }
