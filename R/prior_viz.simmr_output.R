#' Plot the prior distribution for a simmr run
#'
#' This function takes the output from \code{\link{simmr_mcmc}} or
#' \code{\link{simmr_ffvb}} and plots the prior distribution to enable visual
#' inspection. This can be used by itself or together with
#' \code{\link{posterior_predictive}} to visually evaluate the influence of
#' the prior on the posterior distribution.
#'
#' @param simmr_out A run of the simmr model from \code{\link{simmr_mcmc}} or
#' \code{\link{simmr_ffvb}}
#' @param group Which group to run it for (currently only numeric rather than group names)
#' @param plot Whether to create a density plot of the prior or not. The simulated prior values are returned as part of the object
#' @param include_posterior Whether to include the posterior distribution on top of the priors. Defaults to TRUE
#' @param n_sims The number of simulations from the prior distribution
#' @param scales The type of scale from \code{facet_wrap} allowing for \code{fixed}, \code{free}, \code{free_x}, \code{free_y}
#'
#' @returns A list containing \code{plot}: the ggplot object (useful if requires customisation), and \code{sim}: the simulated prior values which can be compared with the posterior densities
#'
#' @export
#'
#' @examples
#' \donttest{
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
#' prior <- prior_viz(simmr_1_out)
#' head(prior$p_prior_sim)
#' summary(prior$p_prior_sim)
#' }
prior_viz <- function(simmr_out,
                      group = 1,
                      plot = TRUE,
                      include_posterior = TRUE,
                      n_sims = 10000,
                      scales = "free") {
  UseMethod("prior_viz")
}
#' @export
prior_viz.simmr_output <- function(simmr_out,
                                   group = 1,
                                   plot = TRUE,
                                   include_posterior = TRUE,
                                   n_sims = 10000,
                                   scales = "free") {
  if (inherits(simmr_out, "simmr_mcmc_object")) {
    # Can't do more than 1 group for now
    assert_int(group, lower = 1, upper = simmr_out$input$n_groups)

    # Get group_name
    n_groups <- simmr_out$input$n_groups
    plot_title_1 <- "Prior distributions"
    plot_title_2 <- "Prior and posterior distributions"
    if (n_groups > 1) {
      group_name <- levels(simmr_out$input$group)[group]
      plot_title_1 <- paste0(plot_title_1, ": ", group_name)
      plot_title_2 <- paste0(plot_title_2, ": ", group_name)
    }

    # Plot and/or output the prior
    mu_f_mean <- simmr_out$output[[1]]$model$data()$mu_f_mean
    sigma_f_sd <- simmr_out$output[[1]]$model$data()$sigma_f_sd
    n_sources <- simmr_out$input$n_sources

    # Now simulate some ps
    p_prior_sim <- matrix(NA, ncol = n_sources, nrow = n_sims)
    for (i in 1:n_sims) {
      f <- rnorm(n_sources, mean = mu_f_mean, sd = sigma_f_sd)
      p_prior_sim[i, ] <- exp(f) / sum(exp(f))
    }
    colnames(p_prior_sim) <- simmr_out$input$source_names
    if (plot) {
      df <- reshape2::melt(p_prior_sim)
      colnames(df) <- c("Num", "Source", "Proportion")
      df$Type <- "Prior"
      out_all <- simmr_out$output[[group]]$BUGSoutput$sims.list$p
      df2 <- reshape2::melt(out_all)
      colnames(df2) <- c("Num", "Source", "Proportion")
      df2$Type <- "Posterior"
      df_all <- rbind(df2, df)
      if (include_posterior) {
        g <- ggplot(
          df_all,
          aes(
            x = Proportion,
            fill = Source,
            linetype = Type
          )
        ) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(aes(y = after_stat(density)), alpha = 0.5) +
          theme_bw() +
          ggtitle(plot_title_2) +
          ylab("Density") +
          facet_wrap("~ Source", scales = scales)
      } else {
        g <- ggplot(
          df,
          aes(
            x = Proportion,
            fill = Source
          )
        ) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(aes(y = after_stat(density)), alpha = 0.5, linetype = 0) +
          theme_bw() +
          ggtitle(plot_title_1) +
          ylab("Density") +
          facet_wrap("~ Source", scales = scales)
      }
      print(g)
    }

    # Return the simulations
    if (exists("g")) {
      invisible(list(plot = g, p_prior_sim = p_prior_sim))
    } else {
      invisible(p_prior_sim)
    }
  } else if (inherits(simmr_out, "simmr_ffvb_object") == TRUE) {
    # Can't do more than 1 group for now
    assert_int(group, lower = 1, upper = simmr_out$input$n_groups)

    # Get group_name
    n_groups <- simmr_out$input$n_groups
    plot_title_1 <- "Prior distributions"
    plot_title_2 <- "Prior and posterior distributions"
    if (n_groups > 1) {
      group_name <- levels(simmr_out$input$group)[group]
      plot_title_1 <- paste0(plot_title_1, ": ", group_name)
      plot_title_2 <- paste0(plot_title_2, ": ", group_name)
    }

    # Plot and/or output the prior
    mu_f_mean <- simmr_out$output[[1]]$model$data$mu_f_mean
    sigma_f_sd <- simmr_out$output[[1]]$model$data$sigma_f_sd
    n_sources <- simmr_out$input$n_sources

    # Now simulate some ps
    p_prior_sim <- matrix(NA, ncol = n_sources, nrow = n_sims)
    for (i in 1:n_sims) {
      f <- rnorm(n_sources, mean = mu_f_mean, sd = sigma_f_sd)
      p_prior_sim[i, ] <- exp(f) / sum(exp(f))
    }
    colnames(p_prior_sim) <- simmr_out$input$source_names
    if (plot) {
      df <- reshape2::melt(p_prior_sim)
      colnames(df) <- c("Num", "Source", "Proportion")
      df$Type <- "Prior"
      out_all <- simmr_out$output[[group]]$BUGSoutput$sims.list$p
      df2 <- reshape2::melt(out_all)
      colnames(df2) <- c("Num", "Source", "Proportion")
      df2$Type <- "Posterior"
      df_all <- rbind(df2, df)
      if (include_posterior) {
        g <- ggplot(
          df_all,
          aes(
            x = Proportion,
            y = after_stat(density),
            fill = Source,
            linetype = Type
          )
        ) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(alpha = 0.5) +
          theme_bw() +
          ggtitle(plot_title_2) +
          ylab("Density") +
          facet_wrap("~ Source", scales = scales)
      } else {
        g <- ggplot(
          df,
          aes(
            x = Proportion,
            y = after_stat(density),
            fill = Source
          )
        ) +
          scale_fill_viridis(discrete = TRUE) +
          geom_density(alpha = 0.5, linetype = 0) +
          theme_bw() +
          ggtitle(plot_title_1) +
          ylab("Density") +
          facet_wrap("~ Source", scales = scales)
      }
      print(g)
    }

    # Return the simulations
    if (exists("g")) {
      invisible(list(plot = g, p_prior_sim = p_prior_sim))
    } else {
      invisible(p_prior_sim)
    }
  }
}
