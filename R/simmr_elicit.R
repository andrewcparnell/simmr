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
#' procedure will not match the target dietary proportion means and standard
#' deviations exactly. If this problem is severe, try increasing the
#' \code{n_sims} value.
#'
#' @param n_sources The number of sources required
#' @param proportion_means The desired prior proportion means. These should sum
#' to 1. Should be a vector of length \code{n_sources}
#' @param proportion_sds The desired prior proportions standard deviations.
#' These have no restricted sum but should be reasonable estimates for a
#' proportion.
#' @param n_sims The number of simulations for which to run the optimisation
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
#' # Simple full run
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # See simmr_mcmc and vignettes for full example run
#' data(simmr_data_1)
#'
#' # Load into simmr and run
#' simmr_1 = simmr_data_1 %>% 
#'   simmr_load %>% 
#'   simmr_mcmc
#' 
#' # Summary
#' simmr_1 %>% summary('quantiles')
#' # A bit vague:
#' deviance    32.271 33.932 35.554 37.780 44.051
#' Source A     0.029  0.113  0.205  0.310  0.499
#' Source B     0.145  0.233  0.283  0.335  0.448
#' Source C     0.217  0.254  0.274  0.295  0.342
#' Source D     0.032  0.123  0.206  0.301  0.462
#' sd[tracer1]  0.171  0.397  0.517  0.669  1.109
#' sd[tracer2]  0.053  0.275  0.407  0.554  0.929
#' 
#' # Now suppose I had prior information that: 
#' # proportion means = 0.5,0.2,0.2,0.1 
#' # proportion sds = 0.08,0.02,0.01,0.02
#' prior = simmr_elicit(4, 
#'                      proportion_means = c(0.5,0.2,0.2,0.1),
#'                      proportion_sds = c(0.08,0.02,0.01,0.02))
#' 
#' # Now re-run the model with these priors
#' simmr_1a = simmr_data_1 %>% simmr_load %>%
#'   simmr_mcmc(prior_control=list(means = prior$mean,
#'                                 sd = prior$sd))
#' 
#' # Look t the results
#' simmr_1a %>% summary('quantiles')
#'               2.5%    25%    50%    75%  97.5%
#' deviance    33.283 36.349 38.774 41.556 48.457
#' Source A     0.396  0.452  0.486  0.517  0.582
#' Source B     0.177  0.206  0.220  0.235  0.263
#' Source C     0.169  0.196  0.208  0.220  0.243
#' Source D     0.053  0.072  0.085  0.098  0.131
#' sd[tracer1]  0.178  0.381  0.500  0.639  1.030
#' sd[tracer2]  0.188  0.445  0.589  0.765  1.276
#' 
#' # Look at the new priors
#' simmr_1a %>% prior_viz
#' 
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
