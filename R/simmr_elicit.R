#' Function to allow informative prior distribution to be included in simmr
#'
#' The main \code{\link{simmr_mcmc}} function allows for a prior distribution
#' to be set for the dietary proportions. The prior distribution is specified
#' by transforming the dietary proportions using the centralised log ratio
#' (CLR). The \code{\link{simmr_elicit}} and \code{\link{simmr_elicit}}
#' functions allows the user to specify
#' prior means and standard deviations for each of the dietary proportions, and
#' then finds CLR-transformed values suitable for input into
#' \code{\link{simmr_mcmc}}.
#'
#' The function takes the desired proportion means and standard deviations,
#' and fits an optimised least squares to the means and standard deviations in
#' turn to produced CLR-transformed estimates for use in
#' \code{\link{simmr_mcmc}}. Using prior information in SIMMs is highly
#' desirable given the restricted nature of the inference. The prior
#' information might come from previous studies, other experiments, or other
#' observations of e.g. animal behaviour.
#'
#' Due to the nature of the restricted space over which the dietary proportions
#' can span, and the fact that this function uses numerical optimisation, the
#' procedure will not match the target dietary proportion meanss and standard
#' deviations exactly. If this problem is severe, try increasing the
#' \code{n_sims} value.
#'
#' @param n_sources The number of sources required
#' @param proportion_means The desired prior proportion means. These should sum
#' to 1. Should be a vector of length \code{n_sources}
#' @param proportion_sds The desired prior proportions standard deviations.
#' These have no restricted sum but should be reasonable estimates for a
#' proportion.
#' @param n_sims The number of simulations for which to run the optimsiation
#' routine.
#' @return A list object with two components \item{mean }{The best estimates of
#' the mean to use in \code{control.prior} in \code{\link{simmr_mcmc}}}
#' \item{sd }{The best estimates of the standard deviations to use in
#' \code{control.prior} in \code{\link{simmr_mcmc}}}
#'
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' 
#' @importFrom stats rnorm optim
#' @importFrom compositions clrInv
#' @importFrom boot logit
#'
#' @examples
#'
#' \dontrun{
# Data set: 10 observations, 2 tracers, 4 sources
#' mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54, 
#'                -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89, 
#'                11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10)
#' colnames(mix) = c('d13C','d15N')
#' s_names=c('Source A','Source B','Source C','Source D')
#' s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, #' nrow=4)
#' s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4)
#' c_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol=2, nrow=4)
#' c_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol=2, nrow=4)
#' conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4)
#' 
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
#' # MCMC run
#' simmr_1_out = simmr_mcmc(simmr_1)
#' 
#' # Summary
#' summary(simmr_1_out,'quantiles')
#' # A bit vague:
#' #           2.5%   25%   50%   75% 97.5%
#' # Source A 0.029 0.115 0.203 0.312 0.498
#' # Source B 0.146 0.232 0.284 0.338 0.453
#' # Source C 0.216 0.255 0.275 0.296 0.342
#' # Source D 0.032 0.123 0.205 0.299 0.465
#' 
#' # Now suppose I had prior information that: 
#' # proportion means = 0.5,0.2,0.2,0.1 
#' # proportion sds = 0.08,0.02,0.01,0.02
#' prior=simmr_elicit(4, c(0.5,0.2,0.2,0.1),c(0.08,0.02,0.01,0.02))
#' 
#' simmr_1a_out = simmr_mcmc(simmr_1,prior_control=list(means=prior$mean,sd=prior$sd#' ))
#' 
#' summary(simmr_1a_out,'quantiles')
#' # Much more precise:
#' #           2.5%   25%   50%   75% 97.5%
#' # Source A 0.441 0.494 0.523 0.553 0.610
#' # Source B 0.144 0.173 0.188 0.204 0.236
#' # Source C 0.160 0.183 0.196 0.207 0.228
#' # Source D 0.060 0.079 0.091 0.105 0.135 
#' }
#' @export simmr_elicit
simmr_elicit <-
  function(n_sources,
           proportion_means = rep(1/n_sources, n_sources),
           proportion_sds = rep(0.1, n_sources),
           n_sims = 1000) {
    # proportion_means must be a vector of length n_sources
    if (length(proportion_means) != n_sources)
      stop('proportion_means must be of same length as the number of sources')
    # proportion_sds must be a vector of length n_sources
    if (length(proportion_sds) != n_sources)
      stop('proportion_sds must be of same length as the number of sources')
    if(any(proportion_sds==0)) 
      stop("No proportion_sds should be 0 as this will mean that food source is not being consumed.")
    
    
    cat('Running elicitation optimisation routine...\n')
    
    # Perform an optimisation to match the standard deviations to a normal distribution
    clr_opt_mean = function(pars) {
      means = c(pars[(1:(n_sources - 1))], -sum(pars[(1:(n_sources - 1))]))
      # Generate n_sims observations from this normal distribution
      f = matrix(
        stats::rnorm(
          n_sims * n_sources,
          mean = means,
          sd = rep(1, n_sources)
        ),
        ncol = n_sources,
        nrow = n_sims,
        byrow = TRUE
      )
      p = as.matrix(compositions::clrInv(f))
      return(sum((
        boot::logit(apply(p, 2, 'mean')) - boot::logit(proportion_means)
      ) ^ 2))
    }
    
    opt_mean = stats::optim(c(rep(0, (n_sources - 1))), clr_opt_mean, method =
                              "SANN")
    if (opt_mean$convergence != 0) {
      warning(
        "Optimisation for means did not converge properly. Please either increase n_sims or adjust proportion_means and proportion_sds"
      )
    } else {
      cat("Mean optimisation successful.\n")
    }
    best_mean = c(opt_mean$par, -sum(opt_mean$par))
    
    options(warn = -1)
    clr_opt_sd = function(pars) {
      sd = pars[1:(n_sources)]
      # Generate n_sims observations from this normal distribution
      f = matrix(
        stats::rnorm(n_sims * n_sources, mean = best_mean, sd = sd),
        ncol = n_sources,
        nrow = n_sims,
        byrow = TRUE
      )
      p = as.matrix(compositions::clrInv(f))
      return(sum((
        boot::logit(apply(p, 2, 'sd')) - boot::logit(proportion_sds)
      ) ^ 2))
    }
    
    opt_sd = stats::optim(rep(1, (n_sources)), clr_opt_sd, method = "SANN")
    if (opt_sd$convergence != 0) {
      warning(
        "Optimisation for stand deviations did not converge properly. Please either increase n_sims or adjust proportion_means and proportion_sds"
      )
    } else {
      cat("Standard deviation optimisation successful.\n")
    }
    best_sd = opt_sd$par
    options(warn = 0)
    
    best_f = matrix(
      stats::rnorm(n_sims * n_sources, mean = best_mean, sd = best_sd),
      ncol = n_sources,
      nrow = n_sims,
      byrow = TRUE
    )
    best_p = as.matrix(compositions::clrInv(best_f))
    
    cat('Best fit estimates provide proportion means of:\n')
    cat(round(apply(best_p, 2, 'mean'), 3))
    cat('\n')
    cat('... and best fit standard deviations of:\n')
    cat(round(apply(best_p, 2, 'sd'), 3))
    cat('\n')
    
    return(list(mean = best_mean, sd = best_sd))
    
}
