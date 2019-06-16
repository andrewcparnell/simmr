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
#' @importFrom dplyr '%>%'
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
#' simmr_1 %>% summary(type='correlations')
#' simmr_1 %>% summary(type='statistics')
#' 
#' # Can save these if required
#' ans = simmr_1 %>% summary(type=c('quantiles','statistics'))
#' 
#' # Plot
#' simmr_1 %>% plot(type='boxplot')
#' simmr_1 %>% plot(type='histogram')
#' simmr_1 %>% plot(type='density')
#' simmr_1 %>% plot(type='matrix')
#' 
#' # Compare two sources
#' simmr_1 %>% compare_sources(source_names=c('Source A','Source D'))
#' 
#' # Compare multiple sources
#' simmr_1 %>% compare_sources
#' 
#' # Look at influence of prior
#' simmr_1 %>% prior_viz
#' 
#' # Check model fit
#' simmr_1 %>% posterior_predictive
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
#' simmr_2 %>% plot(type = 'boxplot')
#' 
#' # Print it
#' simmr_2 %>% print
#' 
#' # Summary
#' simmr_2 %>% summary(type = 'diagnostics')
#' simmr_2 %>% summary(type = 'quantiles')
#' 
#' # Plot
#' simmr_2 %>% plot(type = 'isospace')
#' simmr_2 %>% plot(type = 'histogram')
#' simmr_2 %>% plot(type = 'matrix')
#' 
#' #####################################################################################
#' 
#' # Data set 2: 3 isotopes (d13C, d15N and d34S), 30 observations, 4 sources
#' data(simmr_data_2)
#' 
#' # Run model
#' simmr_3 = data_2 %>% simmr_load %>% simmr_mcmc
#' 
#' # Plot all three tracers
#' simmr_3 %>% plot(tracers = c(1, 2))
#' simmr_3 %>% plot(tracers = c(1, 3))
#' simmr_3 %>% plot(tracers = c(2, 3))
#' # See vignette('simmr') for fancier axis labels
#' 
#' # Summary
#' simmr_3 %>% summary(type='diagnostics')
#' simmr_3 %>% summary(type='quantiles')
#' simmr_3 %>% summary(type='correlations')
#' 
#' # Plot
#' simmr_3 %>% plot(type='boxplot')
#' simmr_3 %>% plot(type='histogram')
#' simmr_3 %>% plot(type='density')
#' simmr_3 %>% plot(type='matrix')
#' 
#' #####################################################################################
#' 
#' # Data set 4 - identified by Fry (2014) as a failing of SIMMs
#' # See the vignette for more interpreation of these data and the output
#' 
#' # The data
#' data(square_data) 
#' 
#' # Load into simmr - note no corrections or concentrations
#' simmr_4 = square_data %>% simmr_load %>% simmr_mcmc
#' 
#' # Look at convergence diagnostics
#' simmr_4 %>% summary(type = 'diagnostics')
#' 
#' # Plot the matrix plot - should look very bad
#' simmr_4 %>% plot(type = 'matrix') # Look at the massive correlations here
#' 
#' # Interesting to plot prior/posterior visualisation of this
#' simmr_4 %>% prior_viz # Not much information in these data
#' 
#' # Check fit
#' simmr_4 %>% posterior_predictive # Still a good fit though
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

# Throw a warning if more than 1 and less than 4 observations
if(simmr_in$n_obs >1 & simmr_in$n_obs< 4) warning("Less than 4 observations. There will likely be problems with parameter estimation.")
  
  
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
