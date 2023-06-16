#' Estimate correction factors from stable isotope data with known dietary
#' proportions
#'
#' This function runs a slightly different version of the main
#' \code{\link{simmr_mcmc}} function with the key difference that it estimates
#' the correction factors (sometimes called trophic enrichment or trophic
#' discrimination factors; TEFs/TDFs) for a given set of dietary proportions.
#'
#' The idea is that this code can be used for feeding studies where an
#' organism is fed a known proportional diet with a view to estimating
#' the correction factors to be used in a later stable isotope mixing
#' model when the organisms are observed in the field.
#'
#' The main argument of the function is an object created from
#' \code{\link{simmr_load}} which contains mixture data on a number of tracers
#' and food source means and standard deviations. Any correction factors
#' included in this object will be ignored. The known dietary proportions should be provided for each individual (i.e. should be a matrix with the same number of rows as \code{mix}). It is advisable to have multiple different dietary proportion values as part of the feeding experimental design
#'
#' The output of the function is a posterior distribution on the correction
#' factors for each food source. Just like the output from
#' \code{\link{simmr_mcmc}}, this should be checked for convergence. Examples
#' are included below to help assist with this check and further plots
#'
#' If, after running \code{\link{simmr_mcmc_tdf}} the convergence diagnostics in
#' \code{\link{summary.simmr_output_tdf}} are not satisfactory, the values of
#' \code{iter}, \code{burn} and \code{thin} in \code{mcmc_control} should be
#' increased by e.g. a factor of 10.
#'
#' @param simmr_in An object created via the function \code{\link{simmr_load}}
#' @param p The known dietary proportions for the feeding study. Dietary proportions should be given per individual (even if they are all identical)
#' @param prior_control A list of values including arguments named \code{means}
#' and \code{sd} which represent the prior means and standard deviations of the
#' correction factors. These can usually be left at their default values unless
#' you wish to include to include prior information on them.
#' @param mcmc_control A list of values including arguments named \code{iter}
#' (number of iterations), \code{burn} (size of burn-in), \code{thin} (amount
#' of thinning), and \code{n.chain} (number of MCMC chains).
#' @return An object of class \code{simmr_tdf} with two named top-level
#' components:
#' \item{input}{The \code{simmr_input} object given to the
#' \code{simmr_mcmc} function}
#' \item{output}{A set of MCMC chains of class
#' \code{mcmc.list} from the coda package. These can be analysed using
#' \code{\link{summary.simmr_output_tdf}}}.
#'
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#'
#' @seealso \code{\link{simmr_load}} for creating objects suitable for this
#' function, \code{\link{simmr_mcmc}} for estimating dietary proportions,
#' \code{\link{plot.simmr_input}} for creating isospace plots,
#' \code{\link{summary.simmr_output_tdf}} for summarising output
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
#'
#' @examples
#' \donttest{
#' ## Example of estimating TDFs for a simple system with known dietary proportions
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # Assume p = c(0.25, 0.25, 0.25, 0.25)
#'
#' # The data
#' data(simmr_data_1)
#' # Load into simmr
#' simmr_tdf <- with(
#'   simmr_data_1,
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
#' # Plot
#' plot(simmr_tdf)
#'
#' # MCMC run
#' simmr_tdf_out <- simmr_mcmc_tdf(simmr_tdf,
#'   p = matrix(
#'     rep(
#'       1 / simmr_tdf$n_sources,
#'       simmr_tdf$n_sources
#'     ),
#'     ncol = simmr_tdf$n_sources,
#'     nrow = simmr_tdf$n_obs,
#'     byrow = TRUE
#'   )
#' )

#' # Summary
#' summary(simmr_tdf_out, type = "diagnostics")
#' summary(simmr_tdf_out, type = "quantiles")
#'
#' # Now put these corrections back into the model and check the
#' # iso-space plots and dietary output
#' simmr_tdf_2 <- with(
#'   simmr_data_1,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = simmr_tdf_out$c_mean_est,
#'     correction_sds = simmr_tdf_out$c_sd_est,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' # Plot with corrections now
#' plot(simmr_tdf_2)
#'
#' simmr_tdf_2_out <- simmr_mcmc(simmr_tdf_2)
#' summary(simmr_tdf_2_out, type = "diagnostics")
#' plot(simmr_tdf_2_out, type = "boxplot")
#' }
#'
#' @export
simmr_mcmc_tdf <- function(simmr_in,
                           p = matrix(
                             rep(
                               1 / simmr_in$n_sources,
                               simmr_in$n_sources
                             ),
                             ncol = simmr_in$n_sources,
                             nrow = simmr_in$n_obs,
                             byrow = TRUE
                           ),
                           prior_control = list(
                             c_mean_est = rep(
                               2,
                               simmr_in$n_tracers
                             ),
                             c_sd_est = rep(
                               2,
                               simmr_in$n_tracers
                             )
                           ),
                           mcmc_control = list(
                             iter = 10000,
                             burn = 1000,
                             thin = 10,
                             n.chain = 4
                           )) {
  UseMethod("simmr_mcmc_tdf")
}
#' @export
simmr_mcmc_tdf.simmr_input <- function(simmr_in,
                                       p = matrix(
                                         rep(
                                           1 / simmr_in$n_sources,
                                           simmr_in$n_sources
                                         ),
                                         ncol = simmr_in$n_sources,
                                         nrow = simmr_in$n_obs,
                                         byrow = TRUE
                                       ),
                                       prior_control = list(
                                         c_mean_est = rep(
                                           2,
                                           simmr_in$n_tracers
                                         ),
                                         c_sd_est = rep(
                                           2,
                                           simmr_in$n_tracers
                                         )
                                       ),
                                       mcmc_control = list(
                                         iter = 10000,
                                         burn = 1000,
                                         thin = 10,
                                         n.chain = 4
                                       )) {
  # Throw warning if n.chain =1
  if (mcmc_control$n.chain == 1) warning("Running only 1 MCMC chain will cause an error in the convergence diagnostics")

  # Throw a warning if less than 4 observations in a group - 1 is ok as it wil do a solo run
  if (min(table(simmr_in$group)) > 1 & min(table(simmr_in$group)) < 4) warning("At least 1 group has less than 4 observations - either put each observation in an individual group or use informative prior information")

  # Set up the model string
  model_string <- "
model {
  # Likelihood
  for (j in 1:J) {
    for (i in 1:N) {
      y[i,j] ~ dnorm(inprod(p[i,]*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p[i,],q[,j]), 1/var_y[i,j])
      var_y[i,j] <- inprod(pow(p[i,]*q[,j],2),pow(s_sd[,j],2)+pow(c_sd[,j],2))/pow(inprod(p[i,],q[,j]),2)
+ pow(sigma[j],2)
    }

  }

  # Prior on sigma
  for(j in 1:J) { sigma[j] ~ dgamma(0.001, sig_upp) }

  # Priors on c
  for (j in 1:J) {
    for (k in 1:K) {
      c_mean[k,j] <- c_mean_j[j]
      c_sd[k,j] <- c_sd_j[j]
    }
    c_mean_j[j] ~ dgamma(c_mean_est[j], 1)
    c_sd_j[j] ~ dgamma(c_sd_est[j], 1)
  }

}
"

  if (simmr_in$n_groups > 1) stop("TDF calculation currently only works for single group data")

  # Loop through all the groups
  curr_rows <- which(simmr_in$group_int == 1)
  curr_mix <- simmr_in$mixtures[curr_rows, , drop = FALSE]

  # Determine if a single observation or not
  if (nrow(curr_mix) == 1) {
    message("Only 1 mixture value, performing a simmr solo run...\n")
    solo <- TRUE
  } else {
    solo <- FALSE
  }

  # Create data object
  data <- with(simmr_in, list(
    y = curr_mix,
    p = p,
    s_mean = source_means,
    s_sd = source_sds,
    N = nrow(curr_mix),
    J = n_tracers,
    q = concentration_means,
    K = n_sources,
    c_sd_est = prior_control$c_sd_est,
    c_mean_est = prior_control$c_mean_est,
    sig_upp = ifelse(solo, 0.001, 1000)
  ))

  # Run in JAGS
  output <- R2jags::jags(
    data = data,
    parameters.to.save = c("c_mean", "c_sd"),
    model.file = textConnection(model_string),
    n.chains = mcmc_control$n.chain,
    n.iter = mcmc_control$iter,
    n.burnin = mcmc_control$burn,
    n.thin = mcmc_control$thin
  )

  output_all <- vector("list")
  output_all$input <- simmr_in
  output_all$output <- output
  output_all$c_mean_est <- output$BUGSoutput$median$c_mean
  output_all$c_sd_est <- output$BUGSoutput$median$c_sd
  class(output_all) <- c("simmr_output_tdf")

  return(output_all)
}
