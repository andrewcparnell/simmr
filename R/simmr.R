#' simmr: A package for fitting stable isotope mixing models via JAGS in R
#' 
#' This package runs a
#' simple Stable Isotope Mixing Model (SIMM) and is meant as a longer term
#' replacement to the previous function SIAR.. These are used e.g. to infer dietary
#' proportions of organisms consuming various food sources from observations on
#' the stable isotope values taken from the organisms' tissue samples. However
#' SIMMs can also be used in other scenarios, such as in sediment mixing or the
#' composition of fatty acids. The main functions are \code{\link{simmr_load}}
#' and \code{\link{simmr_mcmc}}. The help files contain examples of the use of
#' this package. See also the vignette for a longer walkthrough.
#' 
#' An even longer term replacement for properly running SIMMs is MixSIAR, which
#' allows for more detailed random effects and the inclusion of covariates.
#' 
#' @name simmr
#' @docType package
#' @author Andrew Parnell <andrew.parnell@mu.ie>
#' 
#' @references Andrew C. Parnell, Donald L. Phillips, Stuart Bearhop, Brice X.
#' Semmens, Eric J. Ward, Jonathan W. Moore, Andrew L. Jackson, Jonathan Grey,
#' David J. Kelly, and Richard Inger. Bayesian stable isotope mixing models.
#' Environmetrics, 24(6):387–399, 2013.
#' 
#' Andrew C Parnell, Richard Inger, Stuart Bearhop, and Andrew L Jackson.
#' Source partitioning using stable isotopes: coping with too much variation.
#' PLoS ONE, 5(3):5, 2010.
#' @keywords multivariate
#' @examples
#' 
#' \dontrun{
#' # A first example with 2 tracers (isotopes), 10 observations, and 4 food sources
#' 
#' # Add in the data
#' mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54, 
#' -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89, 
#' 11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10)
#' colnames(mix) = c('d13C','d15N')
#' s_names=c('Source A','Source B','Source C','Source D')
#' s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, nrow=4)
#' s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4)
#' c_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol=2, nrow=4)
#' c_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol=2, nrow=4)
#' conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4)
#' 
#' # Load into simmr
#' simmr_in = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means=c_means,
#'                      correction_sds=c_sds,
#'                      concentration_means = conc)
#' 
#' # Plot
#' plot(simmr_in)
#' 
#' # MCMC run
#' simmr_out = simmr_mcmc(simmr_in)
#' 
#' # Check convergence - values should all be close to 1
#' summary(simmr_out, type = 'diagnostics')
#' 
#' # Look at output
#' summary(simmr_out, type = 'statistics')
#' 
#' # Look at influence of priors
#' prior_viz(simmr_out)
#' 
#' # Plot output
#' plot(simmr_out, type = 'histogram')
#' }
NULL



