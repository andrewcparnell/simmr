#' Run a \code{simmr_input} object through the Fixed Form Variational
#' Vayes(FFVB) function
#'
#' This is the main function of simmr. It takes a \code{simmr_input} object
#' created via \code{\link{simmr_load}}, runs it in fixed form
#' Variational Bayes to determine the dietary proportions, and then
#' outputs a \code{simmr_output} object for further analysis and plotting
#' via \code{\link{summary.simmr_output}} and \code{\link{plot.simmr_output}}.
#'

#' @param simmr_in An object created via the function \code{\link{simmr_load}}
#' @param ffvb_control A list of values including arguments named \code{n_output}
#' (number of rows in theta output), \code{mu_0} (prior for mu),
#' and \code{sigma_0} (prior for sigma).
#' @return An object of class \code{simmr_output} with two named top-level
#' components: \item{input }{The \code{simmr_input} object given to the
#' \code{simmr_ffvb} function} \item{output }{A set of outputs produced by
#' the ffvb function. These can be analysed using the
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
# simmr_1 <- with(
#' #  geese_data_day1,
#' #  simmr_load(
#' #    mixtures = mixtures,
#' #    source_names = source_names,
#' #    source_means = source_means,
#' #    source_sds = source_sds,
#' #    correction_means = correction_means,
#' #    correction_sds = correction_sds,
#' #    concentration_means = concentration_means
#' #  )
# )
#'
#' # Plot
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # MCMC run
#' simmr_1_out <- simmr_ffvb(simmr_1)
#'
#' # Print it
#' print(simmr_1_out)
#'
#' # Summary
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
# simmr_2 <- with(
#' #  geese_data_day1,
#' #  simmr_load(
#' #    mixtures = mixtures[1, , drop = FALSE],
#' #    source_names = source_names,
#' #    source_means = source_means,
#' #    source_sds = source_sds,
#' #    correction_means = correction_means,
#' #    correction_sds = correction_sds,
#' #    concentration_means = concentration_means
#' #  )
# )
#'
#' # Plot
#' plot(simmr_2)
#'
#' # FFVB run - automatically detects the single observation
#' simmr_2_out <- simmr_ffvb(simmr_2)
#'
#' # Print it
#' print(simmr_2_out)
#'
#' # Summary
#' summary(simmr_2_out)
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
# simmr_3 <- with(
#' #  simmr_data_2,
#' #  simmr_load(
#' #    mixtures = mixtures,
#' #    source_names = source_names,
#' #    source_means = source_means,
#' #    source_sds = source_sds,
#' #    correction_means = correction_means,
#' #    correction_sds = correction_sds,
#' #    concentration_means = concentration_means
#' #  )
# )
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
#' # FFVB run
#' simmr_3_out <- simmr_ffvb(simmr_3)
#'
#' # Print it
#' print(simmr_3_out)
#'
#' # Summary
#' summary(simmr_3_out)
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
#' # Data set 4 - identified by Fry (2014) as a failing of SIMMs
#' # See the vignette for more interpreation of these data and the output
#'
#' # The data
#' data(square_data)
# simmr_4 <- with(
#' #  square_data,
#' #  simmr_load(
#' #    mixtures = mixtures,
#' #    source_names = source_names,
#' #    source_means = source_means,
#' #    source_sds = source_sds
#' #  )
# )
#'
#' # Get summary
#' print(simmr_4)
#'
#' # Plot
#' plot(simmr_4)
#'
#' # FFVB run
#' simmr_4_out <- simmr_ffvb(simmr_4)
#'
#' # Print it
#' print(simmr_4_out)
#'
#' # Summary
#' summary(simmr_4_out)
#' ans <- summary(simmr_4_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(simmr_4_out)
#' plot(simmr_4_out, type = "boxplot")
#' plot(simmr_4_out, type = "histogram")
#' plot(simmr_4_out, type = "density")
#' plot(simmr_4_out, type = "matrix") # Look at the massive correlations here
#'
#' #####################################################################################
#'
#' # Data set 5 - Multiple groups Geese data from Inger et al 2006
#'
#' # Do this in raw data format - Note that there's quite a few mixtures!
#' data(geese_data)
# simmr_5 <- with(
#' #  geese_data,
#' #  simmr_load(
#' #    mixtures = mixtures,
#' #    source_names = source_names,
#' #    source_means = source_means,
#' #    source_sds = source_sds,
#' #    correction_means = correction_means,
#' #    correction_sds = correction_sds,
#' #    concentration_means = concentration_means,
#' #    group = groups
#' #  )
# )
#'
#' # Plot
#' plot(simmr_5,
#'   xlab = expression(paste(delta^13, "C (\\u2030)", sep = "")),
#'   ylab = expression(paste(delta^15, "N (\\u2030)", sep = "")),
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
simmr_ffvb <- function(simmr_in,
                       ffvb_control = list(
                         n_output = 3600,
                         mu_0 = 0,
                         sigma_0 = 1
                       )) {
  #### make sure this has right file name
  Rcpp::sourceCpp("~/Documents/GitHub/simmrwithffvb_2/src/run_VB.cpp")

  # Throw a warning if less than 4 observations in a group - 1 is ok as it wil do a solo run
  if (min(table(simmr_in$group)) > 1 & min(table(simmr_in$group)) < 4) warning("At least 1 group has less than 4 observations - either put each observation in an individual group or use informative prior information")


  output <- vector("list", length = simmr_in$n_groups)
  names(output) <- levels(simmr_in$group)
  K <- simmr_in$n_sources
  n_tracers <- simmr_in$n_tracers
  n_output <- ffvb_control$n_output
  mu_a <- ffvb_control$mu_0
  sigma_a <- ffvb_control$sigma_0

  lambdares <- matrix(rep(NA, (((K + (K * (K + 1)) / 2)) + n_tracers * 2) * simmr_in$n_groups),
    nrow = (((K + (K * (K + 1)) / 2)) + n_tracers * 2),
    ncol = simmr_in$n_groups
  )
  thetares <- matrix(rep(NA, ((K + n_tracers) * n_output * simmr_in$n_groups)),
    ncol = (K + n_tracers),
    nrow = n_output * simmr_in$n_groups
  )

  mylist <- vector("list", length = simmr_in$n_groups)

  names(mylist) <- levels(simmr_in$group)

  p_fun <- function(x) exp(x) / sum(exp(x))


  # Loop through all the groups
  for (i in 1:simmr_in$n_groups) {
    if (simmr_in$n_groups > 1) cat(paste("\nRunning for group", levels(simmr_in$group)[i], "\n\n"))

    curr_rows <- which(simmr_in$group_int == i)
    curr_mix <- simmr_in$mixtures[curr_rows, , drop = FALSE]

    # Determine if a single observation or not
    if (nrow(curr_mix) == 1) {
      cat("Only 1 mixture value, performing a simmr solo run...\n")
      solo <- TRUE
    } else {
      solo <- FALSE
    }


    n_tracers <- simmr_in$n_tracers
    n_sources <- simmr_in$n_sources
    s_names <- simmr_in$source_names
    K <- simmr_in$n_sources
    S <- 100
    source_means <- simmr_in$source_means
    source_sds <- simmr_in$source_sds
    correction_means <- simmr_in$correction_means
    correction_sds <- simmr_in$correction_sds
    concentration_means <- simmr_in$concentration_means
    y <- curr_mix




    lambdastart <- c(rep(mu_a, K), rep(sigma_a, (((K * (K + 1)) / 2) + n_tracers * 2)))

    lambdares[, i] <- run_VB_cpp(
      lambdastart, K, n_tracers, concentration_means,
      source_means, correction_means, correction_sds,
      source_sds, y
    )

    thetares[(1 + n_output * (i - 1)):(n_output * i), ] <-
      sim_thetacpp(n_output, lambdares[, i], K, n_tracers)



    p <- t(apply(thetares[(1 + n_output * (i - 1)):(n_output * i), 1:K], 1, p_fun))
    sigma <- 1 / sqrt(thetares[(1 + n_output * (i - 1)):(n_output * i), (K + 1):(K + n_tracers)])
    colnames(p) <- simmr_in$source_names
    colnames(sigma) <- colnames(simmr_in$source_sds)

    mylist[[i]] <- list(
      source_names = simmr_in$source_names,
      theta = thetares,
      groupnames = simmr_in$group,
      lambdares = lambdares,
      BUGSoutput = list(
        sims.list = list(
          p = p,
          sigma = sigma
        ),
        sims.matrix = cbind(
          p,
          sigma
        )
      ),
      model = list(data = list(
        mu_f_mean = c(rep(mu_a, K)),
        sigma_f_sd = c(rep(sigma_a, K))
      ))
    )
  }

  output_all <- vector("list")
  output_all$input <- simmr_in
  output_all$output <- mylist

  class(output_all) <- c("simmr_output", "ffvb")

  return(output_all)
}
