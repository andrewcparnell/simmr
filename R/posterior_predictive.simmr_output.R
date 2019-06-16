#' Plot the posterior predictive distribution for a simmr run
#'
#' This function takes the output from \code{\link{simmr_mcmc}} and plots the posterior predictive distribution to enable visualisation of model fit. The simulated posterior predicted values are returned as part of the object and can be saved for external use
#'
#' @param simmr_out A run of the simmr model from \code{\link{simmr_mcmc}}
#' @param prob The probability interval for the posterior predictives. The default is 0.5 (i.e. 50pc intervals)
#' @param plot_ppc Whether to create a bayesplot of the posterior predictive or not. 
#'
#' @importFrom bayesplot ppc_intervals
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # See simmr_mcmc and vignettes for full example run
#' data(simmr_data_1)
#'
#' # Load into simmr (includes iso-space plot)
#' simmr_1 = simmr_data_1 %>% 
#'   simmr_load %>% 
#'   simmr_mcmc
#'
#' # Prior predictive
#' simmr_1 %>% posterior_predictive
#' }
posterior_predictive = function(simmr_out,
                                prob = 0.5,
                                plot_ppc = TRUE) {
  UseMethod('posterior_predictive')
}  
#' @export
posterior_predictive.simmr_output = function(simmr_out,
                                             prob = 0.5,
                                             plot_ppc = TRUE) {
  
# Get the original jags script
model_string_old = simmr_out$output$model$model()
  
# Plug in y_pred
copy_lines = model_string_old[6]
copy_lines = sub("y\\[i","y_pred\\[i",copy_lines)
model_string_new = c(model_string_old[1:6],copy_lines,model_string_old[7:length(model_string_old)])
  
# Re-Run in JAGS
output = R2jags::jags(data=simmr_out$output$model$data(), 
                      parameters.to.save = c("y_pred"),
                      model.file = textConnection(paste0(model_string_new, collapse = '\n')),
                      n.chains = simmr_out$output$BUGSoutput$n.chains,
                      n.iter = simmr_out$output$BUGSoutput$n.iter,
                      n.burnin = simmr_out$output$BUGSoutput$n.burnin,
                      n.thin = simmr_out$output$BUGSoutput$n.thin)

y_post_pred = output$BUGSoutput$sims.list$y_pred

if(plot_ppc) {
  y_rep = y_post_pred
  dim(y_rep) = c(dim(y_post_pred)[1], dim(y_post_pred)[2]*dim(y_post_pred)[3])
  curr_mix = simmr_out$input$mixtures
  g = ppc_intervals(
    y = as.vector(curr_mix),
    yrep = y_rep,
    x = rep(1:nrow(curr_mix), simmr_out$input$n_tracers),
    prob = prob,
    fatten = 1
  )+ ggplot2::ylab("Tracer value") + ggplot2::xlab('Observation') + ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks = 1:simmr_out$input$n_obs)
  print(g)
}
# Return the simulations
invisible(y_post_pred)

}