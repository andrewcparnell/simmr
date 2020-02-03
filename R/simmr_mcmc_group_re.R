#' Run a \code{simmr_input} object with multiple groups through a simmr Markov chain Monte
#' Carlo (MCMC) function with group as a random effect
#' 
#' This function allows multiple groups to be run through simmr with the group
#' variable treated as a random effect. It takes a \code{simmr_input} object
#' created via \code{\link{simmr_load}}, runs an MCMC to determine the dietary
#' proportions, and then outputs a \code{simmr_output_re} object for further
#' analysis and plotting via \code{\link{summary.simmr_output_re}} and
#' \code{\link{plot.simmr_output_re}}.
#' 
#' If, after running \code{\link{simmr_mcmc_group_re}} the convergence diagnostics in
#' \code{\link{summary.simmr_output_re}} are not satisfactory, the values of
#' \code{iter}, \code{burn} and \code{thin} in \code{mcmc_control} should be
#' increased by a factor of 10.
#' 
#' @param simmr_in An object created via the function \code{\link{simmr_load}}.
#' This object must have multiple groups
#' @param source_control A binary matrix which identifies which sources are to be 
#' included in which groups. The matrix should have the same number of rows 
#' as the number of groups, and the same number of columns as the total number of 
#' sources. If NULL (default) will assume all sources contribute to the mixture
#' for all groups
#' @param scale_re The prior scale value of the Cauchy t-distributed 
#' random effect. For detailed prior control of individual groups use the main
#' \code{\link{simmr_mcmc}} function.
#' @param mcmc_control A list of values including arguments named \code{iter}
#' (number of iterations), \code{burn} (size of burn-in), \code{thin} (amount
#' of thinning), and \code{n.chain} (number of MCMC chains).
#' @return An object of class \code{simmr_output_re} with two named top-level
#' components: \item{input }{The \code{simmr_input} object given to the
#' \code{simmr_mcmc} function} \item{output }{A set of MCMC chains of class
#' \code{mcmc.list} from the coda package. These can be analysed using the
#' \code{\link{summary.simmr_output_re}} and \code{\link{plot.simmr_output_re}}
#' functions.}
#' 
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' 
#' @seealso \code{\link{simmr_load}} for creating objects suitable for this
#' function, \code{\link{plot.simmr_input}} for creating isospace plots,
#' \code{\link{summary.simmr_output_re}} for summarising output, and
#' \code{\link{plot.simmr_output_re}} for plotting output. See 
#' \code{\link{simmr_mcmc}} for the main \code{simmr} MCMC function which
#' does not borrow strength between groups and so runs them individually
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
#' # Multiple groups Geese data from Inger et al 2006
#' 
#' # Do this in raw data format - Note that there's quite a few mixtures!
#' mix = matrix(c(10.22, 10.37, 10.44, 10.52, 10.19, 10.45, 9.91, 11.27, 
#'                9.34, 11.68, 12.29, 11.04, 11.46, 11.73, 12.29, 11.79, 11.49, 
#'                11.73, 11.1, 11.36, 12.19, 11.03, 11.21, 10.58, 11.61, 12.16, 
#'                10.7, 11.47, 12.07, 11.75, 11.86, 12.33, 12.36, 11.13, 10.92, 
#'                12.42, 10.95, 12.28, 11.04, 10.76, 10.99, 10.78, 11.07, 10.2, 
#'                11.67, 7.53, 10.65, 10.58, 11.13, 7.73, 10.79, 10.47, 10.82, 
#'                10.41, 11.1, 10.95, 10.76, 10.83, 10.25, 10.52, 9.94, 9.94, 11.61, 
#'                10.65, 10.76, 11.11, 10.2, 11.27, 10.21, 10.88, 11.21, 11.36, 
#'                10.75, 12.38, 11.16, 11.57, 10.79, 11.13, 10.72, 10.99, 10.38, 
#'                10.95, 10.75, 10.75, 11.05, 10.66, 10.61, 10.9, 11.14, 10.33, 
#'                10.83, 10.75, 9.18, 9.03, 9.05, 8.6, 8.29, 10.32, 10.28, 6.47, 
#'                11.36, 10.75, 11.13, 11.37, 10.86, 10.54, 10.39, 10.66, 9.99, 
#'                11.65, 11.02, 10.67, 8.15, 11.12, 10.95, 11.2, 10.76, 11.32, 
#'                10.85, 11.74, 10.46, 10.93, 12.3, 10.67, 11.51, 10.56, 12.51, 
#'                13.51, 11.98, 12.2, 10.48, 12.4, 13, 11.36, 12.08, 12.39, 12.28, 
#'                12.6, 11.3, 11.1, 11.42, 11.49, 12, 13.35, 11.97, 13.35, 12.75, 
#'                12.55, 12.3, 12.51, 12.61, 10.98, 11.82, 12.27, 12.11, 12.11, 
#'                12.89, 12.99, 12.29, 11.89, 12.74, 12.29, 11.89, 10.56, 9.27, 
#'                10.54, 10.97, 10.46, 10.56, 10.86, 10.9, 11.06, 10.76, 10.64, 
#'                10.94, 10.85, 10.45, 11.15, 11.23, 11.16, 10.94, 11.2, 10.71, 
#'                9.55, 8.6, 9.67, 8.17, 9.81, 10.94, 9.49, 9.46, 7.94, 9.77, 8.07, 
#'                8.39, 8.95, 9.83, 8.51, 8.86, 7.93, 8, 8.33, 8, 9.39, 8.01, 7.59, 
#'                8.26, 9.49, 8.23, 9.1, 8.21, 9.59, 9.37, 9.47, 8.6, 8.23, 8.39, 
#'                8.24, 8.34, 8.36, 7.22, 7.13, 10.64, 8.06, 8.22, 8.92, 9.35, 
#'                7.32, 7.66, 8.09, 7.3, 7.33, 7.33, 7.36, 7.49, 8.07, 8.84, 7.93, 
#'                7.94, 8.74, 8.26, 9.63, 8.85, 7.55, 10.05, 8.23, 7.74, 9.12, 
#'                7.33, 7.54, 8.8, -11.36, -11.88, -10.6, -11.25, -11.66, -10.41, 
#'                -10.88, -14.73, -11.52, -15.89, -14.79, -17.64, -16.97, -17.25, 
#'                -14.77, -15.67, -15.34, -15.53, -17.27, -15.63, -15.94, -14.88, 
#'                -15.9, -17.11, -14.93, -16.26, -17.5, -16.37, -15.21, -15.43, 
#'                -16.54, -15, -16.41, -15.09, -18.06, -16.27, -15.08, -14.39, 
#'                -21.45, -22.52, -21.25, -21.84, -22.51, -21.97, -20.23, -21.64, 
#'                -22.49, -21.91, -21.65, -21.37, -22.9, -21.13, -19.33, -20.29, 
#'                -20.56, -20.87, -21.07, -21.69, -21.17, -21.74, -22.69, -21.06, 
#'                -20.42, -21.5, -20.15, -21.99, -22.3, -21.71, -22.48, -21.86, 
#'                -21.68, -20.97, -21.91, -19.05, -22.78, -22.36, -22.46, -21.52, 
#'                -21.84, -21.3, -21.39, -22.1, -21.59, -20.14, -20.67, -20.31, 
#'                -20.07, -21.2, -20.44, -22.06, -22.05, -21.44, -21.93, -22.47, 
#'                -22.27, -22.19, -22.81, -20.48, -22.47, -18.06, -20.72, -20.97, 
#'                -19.11, -18.4, -20.45, -21.2, -19.74, -20.48, -21.48, -17.81, 
#'                -19.77, -22.56, -14.72, -12.21, -12.35, -13.88, -14.43, -14.65, 
#'                -13.9, -14.12, -10.88, -10.44, -15.33, -13.78, -13.98, -15.22, 
#'                -15.25, -15.76, -15.78, -15.49, -13.02, -15.3, -15.55, -14.35, 
#'                -14.99, -14.83, -16.18, -15.01, -12.87, -14.67, -13.84, -14.89, 
#'                -13.33, -15.04, -14.29, -15.62, -13.99, -15.06, -15.06, -15, 
#'                -14.55, -13.32, -14.34, -14.47, -14.31, -14.18, -16.18, -16.25, 
#'                -15.92, -15.35, -14.29, -15.92, -15.35, -20.22, -21.4, -19.97, 
#'                -20.78, -20.61, -20.58, -20.19, -20.71, -20.59, -20.09, -19.37, 
#'                -20.41, -20.84, -20.75, -20.29, -20.89, -19.69, -20.41, -21.24, 
#'                -19.33, -25.87, -25.4, -27.23, -27.52, -24.55, -17.36, -24.7, 
#'                -27.76, -28.92, -25.98, -26.77, -28.76, -27.7, -24.75, -25.47, 
#'                -26.58, -28.94, -29.13, -26.65, -28.04, -27.5, -29.28, -27.85, 
#'                -27.41, -27.57, -29.06, -25.98, -28.21, -25.27, -14.43, -27.4, 
#'                -27.76, -28.45, -27.35, -28.83, -29.39, -28.86, -28.61, -29.27, 
#'                -20.32, -28.21, -26.3, -28.27, -27.75, -28.55, -27.38, -29.13, 
#'                -28.66, -29.02, -26.04, -26.06, -28.52, -28.51, -27.93, -29.07, 
#'                -28.41, -26.42, -27.71, -27.75, -24.28, -28.43, -25.94, -28, 
#'                -28.59, -22.61, -27.34, -27.35, -29.14), ncol=2, nrow=251)
#' colnames(mix) = c('d13C','d15N')
#' s_names = c("Zostera", "Grass", "U.lactuca", "Enteromorpha")
#' s_means = matrix(c(6.49, 4.43, 11.19, 9.82, -11.17, -30.88, -11.17, 
#'                    -14.06), ncol=2, nrow=4)
#' s_sds = matrix(c(1.46, 2.27, 1.11, 0.83, 1.21, 0.64, 1.96, 1.17), ncol=2, nrow=4)
#' c_means = matrix(c(3.54, 3.54, 3.54, 3.54, 1.63, 1.63, 1.63, 1.63), ncol=2, nrow=4)
#' c_sds = matrix(c(0.74, 0.74, 0.74, 0.74, 0.63, 0.63, 0.63, 0.63), ncol=2, nrow=4)
#' conc = matrix(c(0.03, 0.04, 0.02, 0.01, 0.36, 0.4, 0.21, 0.18), ncol=2, nrow=4)
#' grp = paste('day', (c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
#'         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 
#'         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#'         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 
#'         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 
#'         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
#'         7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
#'         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8)))
#' 
#' # Load this in:
#' simmr_in = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means=c_means,
#'                      correction_sds=c_sds,
#'                      concentration_means = conc,
#'                      group=grp)
#' 
#' # Plot
#' plot(simmr_in,group=1:8,xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
#'      ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
#'      title='Isospace plot of Inger et al Geese data')
#' 
#' # Run MCMC for each group
#' simmr_out = simmr_mcmc_group_re(simmr_in)
#' 
#' # Check for convergence
#' summary(simmr_out, type = 'diagnostics')
#' 
#' # Summarise output
#' summary(simmr_out,type='quantiles')
#' 
#' # Plot - only a single group allowed
#' plot(simmr_out,type='boxplot',title='simmr output boxplots')
#' plot(simmr_out,type='histogram',title='simmr output histograms')
#' 
#' # Compare sources within a group
#' compare_sources(simmr_out,source_names=c('Zostera','U.lactuca'),group=2)
#' compare_sources(simmr_out,group=2)
#' 
#' # Compare between groups
#' compare_groups(simmr_out,source='Zostera',groups=1:2)
#' compare_groups(simmr_out,source='Zostera',groups=1:3)
#' compare_groups(simmr_out,source='U.lactuca',groups=c(4:5,7,2))
#' 
#' #############################################################
#' 
#' # Suppose source 1 Zostera was not available in groups 1:4
#' use_sources = cbind(c(rep(0,4),rep(1,5)), 1, 1, 1)
#' simmr_out_noZ = simmr_mcmc_group_re(simmr_in, source_control = use_sources)
#' 
#' # The summary should show no Zostera in groups 1 to 4
#' summary(simmr_out_noZ,type='statistics')
#' plot(simmr_out_noZ, type = 'histogram')
#' 
#' #############################################################
#' 
#' # Suppose some of the mixture isotope values were missing
#' # (This is a short hack to create missing values in the mixtures)
#' Suppose the first 4 on the 2nd isotope were missing
#' simmr_in_miss = simmr_in
#' simmr_in_miss$mixtures[1:4,2] = NA
#' 
#' # Everything should still run fine
#' simmr_out_miss = simmr_mcmc_group_re(simmr_in_miss)
#' summary(simmr_out_miss,type='statistics')
#' plot(simmr_out_miss, type = 'histogram')
#' 
#' }
#' 
#' @export
simmr_mcmc_group_re = function(simmr_in, 
                               source_control = matrix(1,nrow=simmr_in$n_groups,
                                                       ncol=simmr_in$n_sources),
                               scale_re=1, 
                               mcmc_control=list(iter=10000,
                                                 burn=1000,
                                                 thin=10,
                                                 n.chain=4)) {
  UseMethod('simmr_mcmc_group_re') 
}  
#' @export
simmr_mcmc_group_re.simmr_input = function(simmr_in, 
                                           source_control = matrix(1,nrow=simmr_in$n_groups,
                                                                   ncol=simmr_in$n_sources),
                                           scale_re=1, 
                                           mcmc_control=list(iter=10000,
                                                             burn=1000,
                                                             thin=10,
                                                             n.chain=4)) {
  
# Throw warning if n.chain =1
if(mcmc_control$n.chain==1) warning("Running only 1 MCMC chain will cause an error in the convergence diagnostics")

# Stop if there are less than 3 groups  
if(simmr_in$n_groups < 3) stop("This function should only be used when there are at least 3 groups")
    
# Set up the model string
model_string = "
model{
  # Likelihood
  for (j in 1:J) {
    for (i in 1:N) {
      y[i,j] ~ dnorm(inprod(p[group[i],]*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p[group[i],],q[,j]), 1/var_y[i,j])
      var_y[i,j] <- inprod(pow(p[group[i],]*q[,j],2),pow(s_sd[,j],2)+pow(c_sd[,j],2))/pow(inprod(p[group[i],],q[,j]),2)
+ pow(sigma[j],2)
    }
  }

  # Prior on sigma
  for(j in 1:J) { sigma[j] ~ dt(0,100,1)T(0,) }

  # CLR prior on p
  for(l in 1:n_groups) {
    p[l,1:K] <- source_control[l,1:K]*expf[l,]/sum(expf[l,])
    for(k in 1:K) {
      expf[l,k] <- exp(f[l,k])
      f[l,k] ~ dnorm(mu_f_mean[k],sigma_re[k]^-2)
    }
  }
  
  # Create overall group proportions
  p_all[1:K] <- expf_all/sum(expf_all)
  for(k in 1:K) { 
    mu_f_mean[k] ~ dnorm(0, 2^-2)
    sigma_re[k] ~ dt(0, sigma_scale, 1)T(0,)
    expf_all[k] <- exp(mu_f_mean[k]) 
  }
}
"

# Create data object
data = with(simmr_in, list(
  y=simmr_in$mixtures,
  s_mean=source_means,
  s_sd=source_sds,
  N=nrow(simmr_in$mixtures),
  J=n_tracers,
  c_mean=correction_means,
  c_sd=correction_sds,
  q=concentration_means,
  K=n_sources,
  source_control = source_control,
  group = group_int,
  n_groups = n_groups,
  sigma_scale = scale_re
))
  
# Run in JAGS
output = R2jags::jags(data=data, 
                      parameters.to.save = c("p", "sigma", "p_all", "sigma_re"),
                      model.file = textConnection(model_string),
                      n.chains = mcmc_control$n.chain,
                      n.iter = mcmc_control$iter,
                      n.burnin = mcmc_control$burn,
                      n.thin = mcmc_control$thin)
  
# Get the names right and interpretable everywhere
# Set the posterior names right
n_tracers = simmr_in$n_tracers
n_sources = simmr_in$n_sources
n_groups = simmr_in$n_groups
t_names = colnames(simmr_in$mixtures)
s_names = simmr_in$source_names
g_names = levels(simmr_in$group)

# Rename the matrix
count = 2
for(i in 1:n_sources) {
  for(j in 1:n_groups) {
    colnames(output$BUGSoutput$sims.matrix)[count] = paste0('p[',g_names[j],',',s_names[i],']')
    count = count + 1
  }
}
for(i in 1:n_sources) {
  colnames(output$BUGSoutput$sims.matrix)[count] = paste0('p_overall[',s_names[i],']')
  count = count + 1
}
for(i in 1:n_tracers) {
  colnames(output$BUGSoutput$sims.matrix)[count] = paste0('sigma[',t_names[i],']')
  count = count + 1
}
for(i in 1:n_sources) {
  colnames(output$BUGSoutput$sims.matrix)[count] = paste0('sigma_re[',s_names[i],']')
  count = count + 1
}

# And the summary
rownames(output$BUGSoutput$summary) = colnames(output$BUGSoutput$sims.matrix)

output_all = list(
  input = simmr_in,
  output = output)
class(output_all) = 'simmr_output_re'

return(output_all)

}
