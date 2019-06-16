#' simmr: A package for fitting stable isotope mixing models via JAGS in R
#' 
#' This package runs a
#' simple Stable Isotope Mixing Model (SIMM) and is meant as a longer term
#' replacement to the previous function SIAR. These are used to infer dietary
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
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
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
#' 
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
#' # Check convergence
#' simmr_1 %>% summary(type = 'convergence')
#' 
#' # Check the priors
#' simmr_1 %>% prior_viz
#' 
#' # Plot the output
#' simmr_1 %>% plot(type = 'boxplot')
#' 
#' # Check the model fit
#' simmr_1 %>% posterior_predictive
#' 
#' }
NULL



