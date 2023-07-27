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
#' \donttest{
#' # Data set: 10 observations, 2 tracers, 4 sources
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Look at the prior influence
#' prior_viz(simmr_1_out)
#'
#' # Summary
#' summary(simmr_1_out, "quantiles")
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
#' prior <- simmr_elicit(4, c(0.5, 0.2, 0.2, 0.1), c(0.08, 0.02, 0.01, 0.02))
#'
#' simmr_1a_out <- simmr_mcmc(simmr_1, prior_control = 
#' list(means = prior$mean, 
#'       sd = prior$sd, 
#'       shape = 1, 
#'       rate = 1))
#'
#' #' # Look at the prior influence now
#' prior_viz(simmr_1a_out)
#'
#' summary(simmr_1a_out, "quantiles")
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
           proportion_means = rep(1 / n_sources, n_sources),
           proportion_sds = rep(0.1, n_sources),
           n_sims = 1000) {
    # proportion_means must be a vector of length n_sources
    assert_numeric(proportion_means,
      len = n_sources,
      lower = 0, upper = 1
    )
    # proportion_sds must be a vector of length n_sources
    assert_numeric(proportion_sds,
      len = n_sources,
      lower = .Machine$double.eps, upper = 1
    )

    low_cis <- proportion_means - 2 * proportion_sds
    high_cis <- proportion_means + 2 * proportion_sds
    if (any(low_cis < 0) | any(high_cis > 1)) warning("Some proportion sds are large and lie at the edge of the simplex. Check results are reasonable by running prior_viz on the object created by simmr_mcmc.")

    cat("Running elicitation optimisation routine...\n")

    # Perform an optimisation to match the standard deviations to a normal distribution
    clr_opt_mean <- function(pars) {
      means <- c(pars[(1:(n_sources - 1))], -sum(pars[(1:(n_sources - 1))]))
      # Generate n_sims observations from this normal distribution
      suppressWarnings(
      f <- matrix(
        stats::rnorm(
          n_sims * n_sources,
          mean = means,
          sd = rep(1, n_sources)
        ),
        ncol = n_sources,
        nrow = n_sims,
        byrow = TRUE
      ), classes = "warning"
      )
      p <- as.matrix(compositions::clrInv(f))
      return(sum((
        boot::logit(apply(p, 2, "mean")) - boot::logit(proportion_means)
      )^2))
    }

    opt_mean <- stats::optim(c(rep(0, (n_sources - 1))), clr_opt_mean,
      method =
        "SANN"
    )
    if (opt_mean$convergence != 0) {
      warning(
        "Optimisation for means did not converge properly. Please either increase n_sims or adjust proportion_means and proportion_sds"
      )
    } else {
      message("Mean optimisation successful.\n")
    }
    best_mean <- c(opt_mean$par, -sum(opt_mean$par))
suppressWarnings({
    clr_opt_sd <- function(pars) {
      sd <- pars[1:(n_sources)]
      # Generate n_sims observations from this normal distribution
      f <- matrix(
        stats::rnorm(n_sims * n_sources, mean = best_mean, sd = sd),
        ncol = n_sources,
        nrow = n_sims,
        byrow = TRUE
      )
      p <- as.matrix(compositions::clrInv(f))
      return(sum((
        boot::logit(apply(p, 2, "sd")) - boot::logit(proportion_sds)
      )^2))
    }

    opt_sd <- stats::optim(rep(1, (n_sources)), clr_opt_sd, method = "SANN")
    if (opt_sd$convergence != 0) {
      warning(
        "Optimisation for stand deviations did not converge properly. Please either increase n_sims or adjust proportion_means and proportion_sds"
      )
    } else {
      message("Standard deviation optimisation successful.\n")
    }
    best_sd <- opt_sd$par
})

    best_f <- matrix(
      stats::rnorm(n_sims * n_sources, mean = best_mean, sd = best_sd),
      ncol = n_sources,
      nrow = n_sims,
      byrow = TRUE
    )

    best_p <- as.matrix(compositions::clrInv(best_f))

    message("Best fit estimates provide proportion means of:\n")
    message(round(apply(best_p, 2, "mean"), 3))
    message("\n")
    message("... and best fit standard deviations of:\n")
    message(round(apply(best_p, 2, "sd"), 3))
    message("\n")
    message("Check these match the input values before proceeding with a model run.\n")

    return(list(mean = best_mean, sd = best_sd))
  }
