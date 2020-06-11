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
#' post_pred = posterior_predictive(simmr_1_out)
#' }
posterior_predictive = function(simmr_out,
                                group = 1,
                                prob = 0.5,
                                plot_ppc = TRUE) {
  UseMethod('posterior_predictive')
}  
#' @export
posterior_predictive.simmr_output = function(simmr_out,
                                             group = 1,
                                             prob = 0.5,
                                             plot_ppc = TRUE) {
  
# Can't do more than 1 group for now
if(length(group) > 1) stop("Multiple groups not supported")
# Get the original jags script
model_string_old = simmr_out$output[[group]]$model$model()
  
# Plug in y_pred
copy_lines = model_string_old[6]
copy_lines = sub("y\\[i","y_pred\\[i",copy_lines)
model_string_new = c(model_string_old[1:6],copy_lines,model_string_old[7:length(model_string_old)])
  
# Re-Run in JAGS
output = R2jags::jags(data=simmr_out$output[[group]]$model$data(), 
                      parameters.to.save = c("y_pred"),
                      model.file = textConnection(paste0(model_string_new, collapse = '\n')),
                      n.chains = simmr_out$output[[group]]$BUGSoutput$n.chains,
                      n.iter = simmr_out$output[[group]]$BUGSoutput$n.iter,
                      n.burnin = simmr_out$output[[group]]$BUGSoutput$n.burnin,
                      n.thin = simmr_out$output[[group]]$BUGSoutput$n.thin)

y_post_pred = output$BUGSoutput$sims.list$y_pred

# Make is look nicer
low_prob = 0.5 - prob/2
high_prob = 0.5 + prob/2
y_post_pred_ci = apply(y_post_pred, 
                       2:3, 
                       'quantile', 
                       prob = c(low_prob, high_prob))
y_post_pred_out = data.frame(
  interval = matrix(y_post_pred_ci, 
                    ncol = simmr_out$input$n_tracers,
                    byrow = TRUE),
  data = as.vector(simmr_out$input$mixtures[simmr_out$input$group_int == group,])
)

y_post_pred_out$outside = y_post_pred_out[,3] > y_post_pred_out[,2] | 
  y_post_pred_out[,3] < y_post_pred_out[,1]
prop_outside = mean(y_post_pred_out$outside)

if(plot_ppc) {
  y_rep = y_post_pred
  dim(y_rep) = c(dim(y_post_pred)[1], dim(y_post_pred)[2]*dim(y_post_pred)[3])
  curr_rows = which(simmr_out$input$group_int==group)  
  curr_mix = simmr_out$input$mixtures[curr_rows,,drop=FALSE]
  g = ppc_intervals(
    y = as.vector(curr_mix),
    yrep = y_rep,
    x = rep(1:nrow(curr_mix), simmr_out$input$n_tracers),
    prob = prob,
    fatten = 1
  ) + ggplot2::ylab("Tracer value") + 
    ggplot2::xlab('Observation') + 
    ggplot2::ggtitle(paste0(prob*100, '% posterior predictive')) + 
    ggplot2::theme_bw() + 
    ggplot2::scale_x_continuous(breaks = 1:simmr_out$input$n_obs)
  print(g)
}
# Return the simulations
invisible(list(table = y_post_pred_out, 
               prop_outside = prop_outside))

}