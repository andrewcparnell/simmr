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
#' @param prior_control A list of values including arguments named: \code{means}
#' and \code{sd} which represent the prior means and standard deviations of the
#' dietary proportions in centralised log-ratio space; \code{shape} and 
#' \code{rate} which represent the prior distribution on the residual standard
#' deviation. These can usually be
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
#' # Plot
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Print it
#' print(simmr_1_out)
#'
#' # Summary
#' summary(simmr_1_out, type = "diagnostics")
#' summary(simmr_1_out, type = "correlations")
#' summary(simmr_1_out, type = "statistics")
#' ans <- summary(simmr_1_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(simmr_1_out, type = "boxplot")
#' plot(simmr_1_out, type = "histogram")
#' plot(simmr_1_out, type = "density")
#' plot(simmr_1_out, type = "matrix")
#'
#' # Compare two sources
#' compare_sources(simmr_1_out, source_names = c("Zostera", "Enteromorpha"))
#'
#' # Compare multiple sources
#' compare_sources(simmr_1_out)
#'
#' #####################################################################################
#'
#' # A version with just one observation
#' data(geese_data_day1)
#' simmr_2 <- with(
#'   geese_data_day1,
#'   simmr_load(
#'     mixtures = mixtures[1, , drop = FALSE],
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
#' plot(simmr_2)
#'
#' # MCMC run - automatically detects the single observation
#' simmr_2_out <- simmr_mcmc(simmr_2)
#'
#' # Print it
#' print(simmr_2_out)
#'
#' # Summary
#' summary(simmr_2_out)
#' summary(simmr_2_out, type = "diagnostics")
#' ans <- summary(simmr_2_out, type = c("quantiles"))
#'
#' # Plot
#' plot(simmr_2_out)
#' plot(simmr_2_out, type = "boxplot")
#' plot(simmr_2_out, type = "histogram")
#' plot(simmr_2_out, type = "density")
#' plot(simmr_2_out, type = "matrix")
#'
#' #####################################################################################
#'
#' # Data set 2: 3 isotopes (d13C, d15N and d34S), 30 observations, 4 sources
#' data(simmr_data_2)
#' simmr_3 <- with(
#'   simmr_data_2,
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
#' # Get summary
#' print(simmr_3)
#'
#' # Plot 3 times
#' plot(simmr_3)
#' plot(simmr_3, tracers = c(2, 3))
#' plot(simmr_3, tracers = c(1, 3))
#' # See vignette('simmr') for fancier axis labels
#'
#' # MCMC run
#' simmr_3_out <- simmr_mcmc(simmr_3)
#'
#' # Print it
#' print(simmr_3_out)
#'
#' # Summary
#' summary(simmr_3_out)
#' summary(simmr_3_out, type = "diagnostics")
#' summary(simmr_3_out, type = "quantiles")
#' summary(simmr_3_out, type = "correlations")
#'
#' # Plot
#' plot(simmr_3_out)
#' plot(simmr_3_out, type = "boxplot")
#' plot(simmr_3_out, type = "histogram")
#' plot(simmr_3_out, type = "density")
#' plot(simmr_3_out, type = "matrix")
#'
#' #####################################################################################
#'
#' # Data set 5 - Multiple groups Geese data from Inger et al 2006
#'
#' # Do this in raw data format - Note that there's quite a few mixtures!
#' data(geese_data)
#' simmr_5 <- with(
#'   geese_data,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means,
#'     group = groups
#'   )
#' )
#'
#' # Plot
#' plot(simmr_5,
#'   xlab = expression(paste(delta^13, "C (per mille)", sep = "")),
#'   ylab = expression(paste(delta^15, "N (per mille)", sep = "")),
#'   title = "Isospace plot of Inger et al Geese data"
#' )
#'
#' # Run MCMC for each group
#' simmr_5_out <- simmr_mcmc(simmr_5)
#'
#' # Summarise output
#' summary(simmr_5_out, type = "quantiles", group = 1)
#' summary(simmr_5_out, type = "quantiles", group = c(1, 3))
#' summary(simmr_5_out, type = c("quantiles", "statistics"), group = c(1, 3))
#'
#' # Plot - only a single group allowed
#' plot(simmr_5_out, type = "boxplot", group = 2, title = "simmr output group 2")
#' plot(simmr_5_out, type = c("density", "matrix"), grp = 6, title = "simmr output group 6")
#'
#' # Compare sources within a group
#' compare_sources(simmr_5_out, source_names = c("Zostera", "U.lactuca"), group = 2)
#' compare_sources(simmr_5_out, group = 2)
#'
#' # Compare between groups
#' compare_groups(simmr_5_out, source = "Zostera", groups = 1:2)
#' compare_groups(simmr_5_out, source = "Zostera", groups = 1:3)
#' compare_groups(simmr_5_out, source = "U.lactuca", groups = c(4:5, 7, 2))
#' }
#'
#' @export
simmr_mcmc <- function(simmr_in,
                       prior_control = list(
                         means = rep(
                           0,
                           simmr_in$n_sources
                         ),
                         sd = rep(
                           1,
                           simmr_in$n_sources
                         ),
                         sigma_shape = rep(3, simmr_in$n_tracers),
                         sigma_rate = rep(3/50, simmr_in$n_tracers)

                       ),
                       mcmc_control = list(
                         iter = 10000,
                         burn = 1000,
                         thin = 10,
                         n.chain = 4
                       )) {
  UseMethod("simmr_mcmc")
}
#' @export
simmr_mcmc.simmr_input <- function(simmr_in,
                                   prior_control = list(
                                     means = rep(0, simmr_in$n_sources),
                                     sd = rep(1, simmr_in$n_sources),
                                     sigma_shape = rep(3, simmr_in$n_tracers),
                                     sigma_rate = rep(3/50, simmr_in$n_tracers)

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
  jags_file <- system.file("jags_models", "mcmc.jags", package = "simmr")

  output <- vector("list", length = simmr_in$n_groups)
  names(output) <- levels(simmr_in$group)

  # Loop through all the groups
  for (i in 1:simmr_in$n_groups) {
    if (simmr_in$n_groups > 1) message("\nRunning for group ", levels(simmr_in$group)[i], "\n\n")

    curr_rows <- which(simmr_in$group_int == i)
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
      s_mean = source_means,
      s_sd = source_sds,
      N = nrow(curr_mix),
      J = n_tracers,
      c_mean = correction_means,
      c_sd = correction_sds,
      q = concentration_means,
      K = n_sources,
      mu_f_mean = prior_control$means,
      sigma_f_sd = prior_control$sd,
      sigma_shape = prior_control$sigma_shape,
      sigma_rate = prior_control$sigma_rate,
      not_solo = ifelse(solo, 0, 1)
    ))

    # Run in JAGS
    output[[i]] <- R2jags::jags(
      data = data,
      parameters.to.save = c("p", "sigma"),
      model.file = jags_file,
      n.chains = mcmc_control$n.chain,
      n.iter = mcmc_control$iter,
      n.burnin = mcmc_control$burn,
      n.thin = mcmc_control$thin
    )

    # Get the names right and interpretable everywhere
    # Set the posterior names right
    n_tracers <- simmr_in$n_tracers
    n_sources <- simmr_in$n_sources
    s_names <- simmr_in$source_names
    colnames(output[[i]]$BUGSoutput$sims.matrix)[2:(2 + n_sources - 1)] <- s_names
    colnames(output[[i]]$BUGSoutput$sims.list$p) <- s_names
    # Also do it in the summary
    rownames(output[[i]]$BUGSoutput$summary)[2:(2 + n_sources - 1)] <- s_names
    sd_names <- paste0("sd[", colnames(simmr_in$mixtures), "]")
    colnames(output[[i]]$BUGSoutput$sims.matrix)[(n_sources + 2):(n_sources + n_tracers + 1)] <- sd_names
    colnames(output[[i]]$BUGSoutput$sims.list$sigma) <- sd_names
    rownames(output[[i]]$BUGSoutput$summary)[(n_sources + 2):(n_sources + n_tracers + 1)] <- sd_names
  }

  output_all <- vector("list")
  output_all$input <- simmr_in
  output_all$output <- output
  class(output_all) <- c("simmr_output", "simmr_mcmc_object")

  return(output_all)
}
