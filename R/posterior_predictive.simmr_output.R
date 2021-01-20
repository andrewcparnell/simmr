#' Plot the posterior predictive distribution for a simmr run
#'
#' This function takes the output from \code{\link{simmr_mcmc}} and plots the posterior predictive distribution to enable visualisation of model fit. The simulated posterior predicted values are returned as part of the object and can be saved for external use
#'
#' @param simmr_out A run of the simmr model from \code{\link{simmr_mcmc}}
#' @param group Which group to run it for (currently only numeric rather than group names)
#' @param prob The probability interval for the posterior predictives. The default is 0.5 (i.e. 50pc intervals)
#' @param plot_ppc Whether to create a bayesplot of the posterior predictive or not.
#'
#' @importFrom bayesplot ppc_intervals
#'
#' @export
#'
#' @examples
#' \dontrun{
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
#' # Print
#' simmr_1
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Prior predictive
#' post_pred <- posterior_predictive(simmr_1_out)
#' }
posterior_predictive <- function(simmr_out,
                                 group = 1,
                                 prob = 0.5,
                                 plot_ppc = TRUE) {
  UseMethod("posterior_predictive")
}
#' @export
posterior_predictive.simmr_output <- function(simmr_out,
                                              group = 1,
                                              prob = 0.5,
                                              plot_ppc = TRUE) {

  # Can't do more than 1 group for now
  if (length(group) > 1) stop("Multiple groups not supported")
  # Get the original jags script
  model_string_old <- simmr_out$output[[group]]$model$model()

  # Plug in y_pred
  copy_lines <- model_string_old[6]
  copy_lines <- sub("y\\[i", "y_pred\\[i", copy_lines)
  model_string_new <- c(model_string_old[1:6], copy_lines, model_string_old[7:length(model_string_old)])

  # Re-Run in JAGS
  output <- R2jags::jags(
    data = simmr_out$output[[group]]$model$data(),
    parameters.to.save = c("y_pred"),
    model.file = textConnection(paste0(model_string_new, collapse = "\n")),
    n.chains = simmr_out$output[[group]]$BUGSoutput$n.chains,
    n.iter = simmr_out$output[[group]]$BUGSoutput$n.iter,
    n.burnin = simmr_out$output[[group]]$BUGSoutput$n.burnin,
    n.thin = simmr_out$output[[group]]$BUGSoutput$n.thin
  )

  y_post_pred <- output$BUGSoutput$sims.list$y_pred

  # Make is look nicer
  low_prob <- 0.5 - prob / 2
  high_prob <- 0.5 + prob / 2
  y_post_pred_ci <- apply(y_post_pred,
    2:3,
    "quantile",
    prob = c(low_prob, high_prob)
  )
  y_post_pred_out <- data.frame(
    interval = matrix(y_post_pred_ci,
      ncol = simmr_out$input$n_tracers,
      byrow = TRUE
    ),
    data = as.vector(simmr_out$input$mixtures[simmr_out$input$group_int == group, ])
  )

  y_post_pred_out$outside <- y_post_pred_out[, 3] > y_post_pred_out[, 2] |
    y_post_pred_out[, 3] < y_post_pred_out[, 1]
  prop_outside <- mean(y_post_pred_out$outside)

  if (plot_ppc) {
    y_rep <- y_post_pred
    dim(y_rep) <- c(dim(y_post_pred)[1], dim(y_post_pred)[2] * dim(y_post_pred)[3])
    curr_rows <- which(simmr_out$input$group_int == group)
    curr_mix <- simmr_out$input$mixtures[curr_rows, , drop = FALSE]
    g <- ppc_intervals(
      y = as.vector(curr_mix),
      yrep = y_rep,
      x = rep(1:nrow(curr_mix), simmr_out$input$n_tracers),
      prob = prob,
      fatten = 1
    ) + ggplot2::ylab("Tracer value") +
      ggplot2::xlab("Observation") +
      ggplot2::ggtitle(paste0(prob * 100, "% posterior predictive")) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(breaks = 1:simmr_out$input$n_obs)
    print(g)
  }
  # Return the simulations
  invisible(list(
    table = y_post_pred_out,
    prop_outside = prop_outside
  ))
}
