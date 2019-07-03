#' Plot the prior distribution for a simmr run
#'
#' This function takes the output from \code{\link{simmr_mcmc}} and plots the prior distribution to enable visual inspection. This can be used by itself or as part of \code{\link{posterior_predictive}} to visually evaluate the influence of the prior on the posterior distribution.
#' 
#' @param simmr_out A run of the simmr model from \code{\link{simmr_mcmc}}
#' @param plot Whether to create a density plot of the prior or not. The simulated prior values are returned as part of the object
#' @param include_posterior Whether to include the posterior distribution on top of the priors. Defaults to TRUE
#' @param n_sims The number of simulations from the prior distribution
#'
#' @export
#' 
#' @examples 
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # See simmr_mcmc and vignettes for full example run
#' data(simmr_data_1)
#'
#' # Load into simmr and run
#' simmr_1 = simmr_data_1 %>% 
#'   simmr_load %>% 
#'   simmr_mcmc
#'
#' # Prior and posterior distributions
#' prior = simmr_1 %>% prior_viz
#' head(prior)
#' summary(prior)
prior_viz = function(simmr_out,
                     plot = TRUE,
                     include_posterior = TRUE,
                     n_sims = 10000) {
  UseMethod('prior_viz')
}  
#' @export
prior_viz.simmr_output = function(simmr_out,
                                  plot = TRUE,
                                  include_posterior = TRUE,
                                  n_sims = 10000) {


# Plot and/or output the prior  
mu_f_mean = simmr_out$output$model$data()$mu_f_mean
sigma_f_sd = simmr_out$output$model$data()$sigma_f_sd
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
  out_all = simmr_out$output$BUGSoutput$sims.list$p
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
      theme(legend.position = 'bottom') + 
      guides(fill=FALSE) + 
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
      theme(legend.position = 'bottom') + 
      guides(fill=FALSE) + 
      ggtitle("Prior distributions")  +
      ylab("Density") +
      facet_wrap("~ Source")
  }
print(g)
}
  
# Return the simulations
invisible(p_prior_sim)

}