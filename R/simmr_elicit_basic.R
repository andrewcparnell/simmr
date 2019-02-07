#' Function to allow informative prior distribution to be included in simmr
#'
#' The main \code{\link{simmr_mcmc}} function allows for a prior distribution
#' to be set for the dietary proportions. The prior distribution is specified
#' by transforming the dietary proportions using the centralised log ratio
#' (CLR). The \code{\link{simmr_elicit}} and \code{\link{simmr_elicit_basic}}
#' functions allows the user to specify
#' prior means and standard deviations for each of the dietary proportions, and
#' then finds CLR-transformed values suitable for input into
#' \code{\link{simmr_mcmc}}. Th version \code{\link{simmr_elicit_basic}} differs
#' from \code{\link{simmr_elicit}} because it does not require a loaded 
#' \code{simmr} object so can be used for e.g. simulation purposes
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
#' @author Andrew Parnell <andrew.parnell@@ucd.ie>
#' @seealso See \code{\link{simmr_elicit}} to use with an existing \code{simmr}
#' object.
#'
#' @importFrom stats rnorm optim
#' @importFrom compositions clrInv
#' @importFrom boot logit
#'
#' @examples
#'
#' \dontrun{
#' # Suppose we had prior information that:
#' # proportion means = 0.5,0.2,0.2,0.1
#' # proportion sds = 0.08,0.02,0.01,0.02
#' prior=simmr_elicit_basic(4,c(0.5,0.2,0.2,0.1),c(0.08,0.02,0.01,0.02))
#' }
#' @export simmr_elicit
simmr_elicit_basic <-
  function(n_sources,
           proportion_means = rep(1/n_sources, n_sources),
           proportion_sds = rep(0.1, n_sources),
           n_sims = 1000) {
    # proportion_means must be a vector of length n_sources
    if (length(proportion_means) != n_sources)
      stop('proportions_means must be of same length as the number of sources')
    # proportion_sds must be a vector of length n_sources
    if (length(proportion_sds) != n_sources)
      stop('proportions_sds must be of same length as the number of sources')
    
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
