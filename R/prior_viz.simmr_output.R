#' Plot the prior distribution for a simmr run
#'
#' This function takes the output from \code{\link{simmr_mcmc}} and plots the prior distribution to enable visual inspection. This can be used by itself or as part of \code{\link{posterior_predictive}} to visually evaluate the influence of the prior on the posterior distribution.
#' 
#' @param simmr_out A run of the simmr model from \code{\link{simmr_mcmc}}
#' @param group Which group to run it for (currently only numeric rather than group names)
#' @param plot Whether to create a density plot of the prior or not. The simulated prior values are returned as part of the object
#' @param include_posterior Whether to include the posterior distribution on top of the priors. Defaults to TRUE
#' @param n_sims The number of simulations from the prior distribution
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54, 
#' -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89, 
#' 11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10)
#' colnames(mix) = c('d13C','d15N')
#' s_names=c('Source A','Source B','Source C','Source D')
#' s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, nrow=4)
#' s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4)
#' c_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol=2, nrow=4)
#' c_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol=2, nrow=4)
#' conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4)
#' 
#' # Load into simmr
#' simmr_1 = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means=c_means,
#'                      correction_sds=c_sds,
#'                      concentration_means = conc)
#' 
#' # Plot
#' plot(simmr_1)
#' 
#' # Print
#' simmr_1
#' 
#' # MCMC run
#' simmr_1_out = simmr_mcmc(simmr_1)
#'
#' # Prior predictive
#' prior = prior_viz(simmr_1_out)
#' head(prior)
#' summary(prior)
#' }
prior_viz = function(simmr_out,
                     group = 1,
                     plot = TRUE,
                     include_posterior = TRUE,
                     n_sims = 10000) {
  UseMethod('prior_viz')
}  
#' @export
prior_viz.simmr_output = function(simmr_out,
                                  group = 1,
                                  plot = TRUE,
                                  include_posterior = TRUE,
                                  n_sims = 10000) {

  # Can't do more than 1 group for now
  if(length(group) > 1) stop("Multiple groups not supported")
  
  # Get group_name
  n_groups = simmr_out$input$n_groups
  plot_title = "Prior and posterior distributions"
  if(n_groups > 1) {
    group_name = levels(simmr_out$input$group)[group]
    plot_title = paste0(plot_title, ': ', group_name)
  }
  
  # Plot and/or output the prior  
  mu_f_mean = simmr_out$output[[1]]$model$data()$mu_f_mean
  sigma_f_sd = simmr_out$output[[1]]$model$data()$sigma_f_sd
  n_sources = simmr_out$input$n_sources
  
  # Now simulate some ps
  p_prior_sim = matrix(NA, ncol = n_sources, nrow = n_sims) 
  for(i in 1:n_sims) {
    f = rnorm(n_sources,mean = mu_f_mean, sd = sigma_f_sd)
    p_prior_sim[i,] = exp(f)/sum(exp(f))
  }
  colnames(p_prior_sim) = simmr_out$input$source_names
  if(plot) {
    df = reshape2::melt(p_prior_sim)
    colnames(df) = c('Num','Source','Proportion')
    df$Type = "Prior"
    out_all = simmr_out$output[[group]]$BUGSoutput$sims.list$p
    df2 = reshape2::melt(out_all)
    colnames(df2) = c('Num','Source','Proportion')
    df2$Type = "Posterior"
    df_all = rbind(df2, df)
    
    if(include_posterior) {
      g=ggplot(df_all,
               aes_string(x="Proportion",
                          y="..density..",
                          fill="Source",
                          linetype="Type")) +
        scale_fill_viridis(discrete=TRUE) +
        geom_density(alpha = 0.5) +
        theme_bw() +
        ggtitle("Prior and posterior distributions")  +
        ylab("Density") +
        facet_wrap("~ Source")
    } else {
      g=ggplot(df,
               aes_string(x="Proportion",
                          y="..density..",
                          fill="Source")) +
        scale_fill_viridis(discrete=TRUE) +
        geom_density(alpha = 0.5, linetype = 0) +
        theme_bw() +
        ggtitle("Prior distributions")  +
        ylab("Density") +
        facet_wrap("~ Source")
    }
    print(g)
  }
  
  # Return the simulations
  invisible(p_prior_sim)
  
}