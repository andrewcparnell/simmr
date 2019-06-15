#' Run a \code{simmr_input} object through the main simmr Markov chain Monte
#' Carlo (MCMC) function
#' 
#' This is the main function of simmr. It takes a \code{simmr_input} object
#' created via \code{\link{simmr_load}}, runs an MCMC to determine the dietary
#' proportions, and then outputs a \code{simmr_output} object for further
#' analysis and plotting via \code{\link{summary.simmr_output}} and
#' \code{\link{plot.simmr_output}}.
#' 
#' If, after running \code{\link{simmr_mcmc}} the convergence diagnostics in
#' \code{\link{summary.simmr_output}} are not satisfactory, the values of
#' \code{iter}, \code{burn} and \code{thin} in \code{mcmc_control} should be
#' increased by a factor of 10.
#' 
#' @param simmr_in An object created via the function \code{\link{simmr_load}}
#' @param prior_control A list of values including arguments named \code{means}
#' and \code{sd} which represent the prior means and standard deviations of the
#' dietary proportions in centralised log-ratio space. These can usually be
#' left at their default values unless you wish to include to include prior
#' information, in which case you should use the function
#' \code{\link{simmr_elicit}}.
#' @param mcmc_control A list of values including arguments named \code{iter}
#' (number of iterations), \code{burn} (size of burn-in), \code{thin} (amount
#' of thinning), and \code{n.chain} (number of MCMC chains).
#' @return An object of class \code{simmr_output} with two named top-level
#' components: \item{input }{The \code{simmr_input} object given to the
#' \code{simmr_mcmc} function} \item{output }{A set of MCMC chains of class
#' \code{mcmc.list} from the coda package. These can be analysed using the
#' \code{\link{summary.simmr_output}} and \code{\link{plot.simmr_output}}
#' functions.}
#' 
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' 
#' @seealso \code{\link{simmr_load}} for creating objects suitable for this
#' function, \code{\link{plot.simmr_input}} for creating isospace plots,
#' \code{\link{summary.simmr_output}} for summarising output, and
#' \code{\link{plot.simmr_output}} for plotting output.
#' 
#' @references Andrew C. Parnell, Donald L. Phillips, Stuart Bearhop, Brice X.
#' Semmens, Eric J. Ward, Jonathan W. Moore, Andrew L. Jackson, Jonathan Grey,
#' David J. Kelly, and Richard Inger. Bayesian stable isotope mixing models.
#' Environmetrics, 24(6):387–399, 2013.
#' 
#' Andrew C Parnell, Richard Inger, Stuart Bearhop, and Andrew L Jackson.
#' Source partitioning using stable isotopes: coping with too much variation.
#' PLoS ONE, 5(3):5, 2010.
#' 
#' @importFrom R2jags jags
#' 
#' @examples
#' \dontrun{
#' ## See the package vignette for a detailed run through of these 4 examples
#' 
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' 
#' # The data
#' data(simmr_data_1)
#' 
#' # Load into simmr
#' simmr_1 = simmr_data_1 %>% simmr_load %>% simmr_mcmc
#' 
#' # Print
#' simmr_1
#' 
#' # Check convergence diagnostics
#' simmr_1 %>% summary(type = 'diagnostics')
#' 
#' # Other summaries
#' summary(simmr_1,type='correlations')
#' summary(simmr_1,type='statistics')
#' 
#' # Can save these if required
#' ans = summary(simmr_1,type=c('quantiles','statistics'))
#' 
#' # Plot
#' plot(simmr_1,type='boxplot')
#' plot(simmr_1,type='histogram')
#' plot(simmr_1,type='density')
#' plot(simmr_1,type='matrix')
#' 
#' # Compare two sources
#' compare_sources(simmr_1,source_names=c('Source A','Source D'))
#' 
#' # Compare multiple sources
#' compare_sources(simmr_1)
#' 
#' # Look at influence of prior
#' prior_viz(simmr_1)
#' 
#' # Check model fit
#' posterior_predictive(simmr_1)
#' 
#' #####################################################################################
#' 
#' # A version with just one observation
#' data(simmr_data_1)
#' simmr_data_1_single = simmr_data_1
#' simmr_data_1_single$mixtures = simmr_data_1 %>% 
#'   '[['('mixtures') %>% 
#'   '['(1,,drop = FALSE)
#' # drop required to keep it as a matrix
#' 
#' simmr_2 = simmr_data_1_single %>% 
#'   simmr_load %>% 
#'   simmr_mcmc
#' 
#' # Plot
#' plot(simmr_2, type = 'boxplot')
#' 
#' # Print it
#' print(simmr_2)
#' 
#' # Summary
#' summary(simmr_2, type = 'diagnostics')
#' summary(simmr_2, type = 'quantiles')
#' 
#' # Plot
#' plot(simmr_2, type = 'isospace')
#' plot(simmr_2, type = 'histogram')
#' plot(simmr_2, type = 'matrix')
#' 
#' #####################################################################################
#' 
#' # Data set 2: 3 isotopes (d13C, d15N and d34S), 30 observations, 4 sources
#' 
#' # The data
#' #data(data2)
#' data_2 = list(
#' mixtures = matrix(c(-11.67, -12.55, -13.18, -12.6, -11.77, -11.21, -11.45, 
#'                -12.73, -12.49, -10.6, -12.26, -12.48, -13.07, -12.67, -12.26, 
#'                -13.12, -10.83, -13.2, -12.24, -12.85, -11.65, -11.84, -13.26, 
#'                -12.56, -12.97, -12.18, -12.76, -11.53, -12.87, -12.49, 7.79, 
#'                7.85, 8.25, 9.06, 9.13, 8.56, 8.03, 7.74, 8.16, 8.43, 7.9, 8.32, 
#'                7.85, 8.14, 8.74, 9.17, 7.33, 8.06, 8.06, 8.03, 8.16, 7.24, 7.24, 
#'                8, 8.57, 7.98, 7.2, 8.13, 7.78, 8.21, 11.31, 10.92, 11.3, 11, 
#'                12.21, 11.52, 11.05, 11.05, 11.56, 11.78, 12.3, 10.87, 10.35, 
#'                11.66, 11.46, 11.55, 11.41, 12.01, 11.97, 11.5, 11.18, 11.49, 
#'                11.8, 11.63, 10.99, 12, 10.63, 11.27, 11.81, 12.25), ncol=3, nrow=30),
#' source_names = c('Source A', 'Source B', 'Source C', 'Source D'),
#' tracer_names = c('d13C','d15N','d34S'),
#' source_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96, 
#'                    10.35, 7.51, 10.31, 9), ncol=3, nrow=4),
#' source_sds = matrix(c(0.46, 0.39, 0.42, 0.48, 0.44, 0.37, 0.49, 0.47, 0.49, 
#'                  0.42, 0.41, 0.42), ncol=3, nrow=4),
#' correction_means = matrix(c(1.3, 1.58, 0.81, 1.7, 1.73, 1.83, 1.69, 3.2, 0.67, 
#'                    2.99, 3.38, 1.31), ncol=3, nrow=4),
#' correction_sds = matrix(c(0.32, 0.64, 0.58, 0.46, 0.61, 0.55, 0.47, 0.45, 0.34, 
#'                  0.45, 0.37, 0.49), ncol=3, nrow=4),
#' concentration_means = matrix(c(0.05, 0.1, 0.06, 0.07, 0.07, 0.03, 0.07, 0.05, 0.1, 
#'                 0.05, 0.12, 0.11), ncol=3, nrow=4)
#' )
#' 
#' # Run model
#' simmr_3 = data_2 %>% simmr_load
#' 
#' # Plot all three tracers
#' plot(simmr_3, tracers = c(1, 2))
#' plot(simmr_3, tracers = c(1, 3))
#' plot(simmr_3, tracers = c(2, 3))
#' # See vignette('simmr') for fancier axis labels
#' 
#' # MCMC run
#' simmr_3_out = simmr_mcmc(simmr_3)
#' 
#' # Print it
#' print(simmr_3_out)
#' 
#' # Summary
#' summary(simmr_3_out,type='diagnostics')
#' summary(simmr_3_out,type='quantiles')
#' summary(simmr_3_out,type='correlations')
#' 
#' # Plot
#' plot(simmr_3_out,type='boxplot')
#' plot(simmr_3_out,type='histogram')
#' plot(simmr_3_out,type='density')
#' plot(simmr_3_out,type='matrix')
#' 
#' #####################################################################################
#' 
#' # Data set 4 - identified by Fry (2014) as a failing of SIMMs
#' # See the vignette for more interpreation of these data and the output
#' 
#' # The data
#' # data(data_4)
#' data_4 = list(
#' mixtures = matrix(c(-14.65, -16.39, -14.5, -15.33, -15.76, -15.15, -15.73, 
#'                -15.52, -15.44, -16.19, 8.45, 8.08, 7.39, 8.68, 8.23, 7.84, 8.48, 
#'                8.47, 8.44, 8.37), ncol=2, nrow=10),
#' source_names = c('Source A', 'Source B', 'Source C', 'Source D'),
#' source_means = matrix(c(-25, -25, -5, -5, 4, 12, 12, 4), ncol=2, nrow=4),
#' source_sds = matrix(c(1, 1, 1, 1, 1, 1, 1, 1), ncol=2, nrow=4))
#' 
#' # Load into simmr - note no corrections or concentrations
#' simmr_4 = data_4 %>% simmr_load %>% simmr_mcmc
#' 
#' # Look at convergence diagnostics
#' summary(simmr_4, type = 'diagnostics')
#' 
#' # Plot the matrix plot - should look very bad
#' plot(simmr_4, type = 'matrix') # Look at the massive correlations here
#' 
#' # Interesting to plot prior/posterior visualisation of this
#' prior_viz(simmr_4) # Not much information in these data
#' 
#' # Check fit
#' posterior_predictive(simmr_4) # Still a good fit though
#' 
#' }
#' 
#' @export
simmr_mcmc = function(simmr_in, 
                      prior_control=list(means=rep(0,
                                                   simmr_in$n_sources),
                                         sd=rep(1,
                                                simmr_in$n_sources)), 
                      mcmc_control=list(iter=10000,
                                        burn=2000,
                                        thin=8,
                                        n.chain=4)) {
  UseMethod('simmr_mcmc') 
}  
#' @export
simmr_mcmc.simmr_input = function(simmr_in, 
                      prior_control=list(means=rep(0,simmr_in$n_sources),
                                         sd=rep(1,simmr_in$n_sources)), 
                      mcmc_control=list(iter=10000,
                                        burn=2000,
                                        thin=8,
                                        n.chain=4)) {

# Throw warning if n.chain =1
if(mcmc_control$n.chain==1) warning("Running only 1 MCMC chain will cause an error in the convergence diagnostics")

# Throw a warning if less than 4 observations
if(simmr_in$n_obs <4) warning("Less than 4 observations. There will likely be problems with parameter estimation.")
  
  
# Set up the model string
model_string = "
model{
  # Likelihood
  for (j in 1:J) {
    for (i in 1:N) {
      y[i,j] ~ dnorm(inprod(p*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p,q[,j]), 1/var_y[j])
    }
    var_y[j] <- inprod(pow(p*q[,j],2),pow(s_sd[,j],2)+pow(c_sd[,j],2))/pow(inprod(p,q[,j]),2)
+ pow(sigma[j],2)
  }

  # Prior on sigma
  for(j in 1:J) { sigma[j] ~ dunif(0,sig_upp) }

  # CLR prior on p
  p[1:K] <- expf/sum(expf)
  for(k in 1:K) {
    expf[k] <- exp(f[k])
    f[k] ~ dnorm(mu_f_mean[k],1/pow(sigma_f_sd[k],2))
  }
}
"

# Determine if a single observation or not
if(nrow(simmr_in$mixtures)==1) {
  cat('Only 1 mixture value, performing a simmr solo run...\n')
  solo=TRUE
} else {
  solo=FALSE
}

# Create data object
data = with(simmr_in,list(
  y=mixtures,
  s_mean=source_means,
  s_sd=source_sds,
  N=nrow(mixtures),
  J=n_tracers,
  c_mean=correction_means,
  c_sd = correction_sds,
  q=concentration_means,
  K=n_sources,
  mu_f_mean=prior_control$means,
  sigma_f_sd=prior_control$sd,
  sig_upp=ifelse(solo,0.001,1000)))

# Run in JAGS
output = R2jags::jags(data=data, 
                      parameters.to.save = c("p", "sigma"),
                      model.file = textConnection(model_string),
                      n.chains = mcmc_control$n.chain,
                      n.iter = mcmc_control$iter,
                      n.burnin = mcmc_control$burn,
                      n.thin = mcmc_control$thin)

# Get the names right and interpretable everywhere
# Set the posterior names right
n_tracers = simmr_in$n_tracers
n_sources = simmr_in$n_sources
s_names = simmr_in$source_names
colnames(output$BUGSoutput$sims.matrix)[2:(2 + n_sources - 1)] = s_names
colnames(output$BUGSoutput$sims.list$p) = s_names
# Also do it in the summary
rownames(output$BUGSoutput$summary)[2:(2 + n_sources - 1)] = s_names
sd_names = paste0('sd[',colnames(simmr_in$mixtures),']')
colnames(output$BUGSoutput$sims.matrix)[(n_sources + 2):(n_sources + n_tracers + 1)] = sd_names
colnames(output$BUGSoutput$sims.list$sigma) = sd_names
rownames(output$BUGSoutput$summary)[(n_sources + 2):(n_sources + n_tracers + 1)] = sd_names

output_all = vector('list')
output_all$input = simmr_in
output_all$output = output
class(output_all) = 'simmr_output'

return(output_all)

}
