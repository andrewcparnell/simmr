#' Function to load in simmr data and check for errors
#'
#' This function takes in the mixture data, food source means and standard
#' deviations, and (optionally) correction factor means and standard
#' deviations, and concentration proportions. It performs some (non-exhaustive)
#' checking of the data to make sure it will run through simmr. It outputs an
#' object of class \code{simmr_input}.
#'
#' For standard stable isotope mixture modelling, the mixture matrix will
#' contain a row for each individual and a column for each isotopic value.
#' \code{simmr} will allow for any number of isotopes and any number of
#' observations, within computational limits. The source means/sds should be
#' provided for each food source on each isotope. The correction means (usually
#' trophic enrichment factors) can be set as zero if required, and should be of
#' the same shape as the source values. The concentration dependence means
#' should be estimated values of the proportion of each element in the food
#' source in question and should be given in proportion format between 0 and 1.
#' At present there is no means to include concentration standard deviations.
#'
#' @param mixtures The mixture data given as a matrix where the number of rows
#' is the number of observations and the number of columns is the number of
#' tracers (usually isotopes)
#' @param source_names The names of the sources given as a character string
#' @param source_means The means of the source values, given as a matrix where
#' the number of rows is the number of sources and the number of columns is the
#' number of tracers
#' @param source_sds The standard deviations of the source values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers
#' @param correction_means The means of the correction values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. If not provided these are set to 0.
#' @param correction_sds The standard deviations of the correction values,
#' given as a matrix where the number of rows is the number of sources and the
#' number of columns is the number of tracers. If not provided these are set to
#' 0.
#' @param concentration_means The means of the concentration values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. These should be between 0 and 1. If not
#' provided these are all set to 1.
#' @param group A grouping variable. These can be a character or factor variable
#' 
#' @return An object of class \code{simmr_input} with the following elements:
#' \item{mixtures }{The mixture data} \item{source_neams }{Source means}
#' \item{sources_sds }{Source standard deviations} \item{correction_means
#' }{Correction means} \item{correction_sds }{Correction standard deviations}
#' \item{concentration_means }{Concentration dependence means} \item{n_obs
#' }{The number of observations} \item{n_tracers }{The number of
#' tracers/isotopes} \item{n_sources }{The number of sources} \item{n_groups
#' }{The number of groups}
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} for complete examples.
#' @examples
#'
#' # A simple example with 10 observations, 2 tracers and 4 sources
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
#' # Load in with simmr_load
#' simmr_1 = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means=c_means,
#'                      correction_sds=c_sds,
#'                      concentration_means = conc)
#'
#' print(simmr_1)
#'
#' @export simmr_load
simmr_load = function(mixtures,
                      source_names,
                      source_means,
                      source_sds,
                      correction_means = NULL,
                      correction_sds = NULL,
                      concentration_means = NULL,
                      group = NULL) {
  # Function to load in data for simmr and check whether it's appropriate for running through simmr_mcmc
  
  # Go through each object and check that it matches the requirements
  
  # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
  if (!is.matrix(mixtures))
    stop("mixtures object must be a matrix")
  n_obs = nrow(mixtures)
  n_tracers = ncol(mixtures)
  
  # Add column names if they're not there
  if (is.null(colnames(mixtures)))
    colnames(mixtures) = paste0('tracer', 1:n_tracers)
  
  # source_names must be a character vector - the length of it is the number of sources
  if (!is.vector(source_names))
    stop("source_names must be a vector")
  source_names = as.character(source_names)
  n_sources = length(source_names)
  
  # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
  if (length(dim(source_means)) != 2)
    stop("source_means must have two dimensions")
  if (length(dim(source_sds)) != 2)
    stop("source_sds must have two dimensions")
  if (nrow(source_means) != n_sources)
    stop('Number of rows in source_means does not match length(source_names)')
  if (ncol(source_means) != n_tracers)
    stop('Number of columns in source_means does not match ncol(mixtures)')
  if (nrow(source_sds) != n_sources)
    stop('Number of rows in source_sds does not match length(source_names)')
  if (ncol(source_sds) != n_tracers)
    stop('Number of columns in source_sds does not match ncol(mixtures)')
  
  # Check that either neither or both of corrections_means and sds are given
  if (is.null(correction_means) &
      !is.null(correction_sds))
    stop("Both correction_means and correction_sds must be supplied")
  if (!is.null(correction_means) &
      is.null(correction_sds))
    stop("Both correction_means and correction_sds must be supplied")
  
  # correction_means and correction_sds must be matrix and of same dimension as n_sources
  if (!is.null(correction_means)) {
    if (length(dim(correction_means)) != 2)
      stop("correction_means must have two dimensions")
    if (length(dim(correction_sds)) != 2)
      stop("correction_sds must have two dimensions")
    if (nrow(correction_means) != n_sources)
      stop('Number of rows in correction_means does not match length(source_names)')
    if (ncol(correction_means) != n_tracers)
      stop('Number of columns in correction_means does not match ncol(mixtures)')
    if (nrow(correction_sds) != n_sources)
      stop('Number of rows in correction_sds does not match length(source_names)')
    if (ncol(correction_sds) != n_tracers)
      stop('Number of columns in correction_sds does not match ncol(mixtures)')
  } else {
    correction_means = matrix(0, ncol = n_tracers, nrow = n_sources)
    correction_sds = matrix(0, ncol = n_tracers, nrow = n_sources)
  }
  
  # concentration_means must be a matrix where all elements are less than 1
  if (!is.null(concentration_means)) {
    if (length(dim(concentration_means)) != 2)
      stop("concentration_means must have two dimensions")
    if (nrow(concentration_means) != n_sources)
      stop('Number of rows in concentration_means does not match length(source_names)')
    if (ncol(concentration_means) != n_tracers)
      stop('Number of columns in concentration_means does not match ncol(mixtures)')
    if (any(concentration_means > 1))
      stop("Invalid values in concentration_means; must all be less than 1")
    if (any(concentration_means == 0))
      stop("Invalid values in concentration_means; must all be non-zero")
    
  } else {
    concentration_means = matrix(1, ncol = n_tracers, nrow = n_sources)
  }
  
  # Check the groups are the right length and structure if given
  if (!is.null(group)) {
    # Allow for the group to be a factor variable
    
    #if(!is.integer(group)) stop("group variable needs to be of integer type. Perhaps use as.integer?")
    group = factor(group, 
                   levels = unique(group),
                   ordered = TRUE)
    group_int = as.integer(group)
    if (min(group_int) != 1)
      stop("Group integers must start at 1")
    if (!all(diff(sort(unique(group_int))) == 1))
      stop("Group integers must proceed sequentially from 1, 2, 3, etc")
    
    if (length(group) != n_obs)
      stop("Number of group values not equal to number of observations")
  } else {
    group = factor(rep(1, n_obs))
    group_int = as.integer(group)
  }
  n_groups = length(unique(group))
  
  # Prepare output and give class
  out = list(
    mixtures = mixtures,
    source_names = source_names,
    source_means = source_means,
    source_sds = source_sds,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means,
    group = group,
    group_int = group_int,
    n_obs = n_obs,
    n_tracers = n_tracers,
    n_sources = n_sources,
    n_groups = n_groups
  )
  
  # Look through to see whether there are any missing values in anything
  if (any(unlist(lapply(out, 'is.na')))) {
    warning("Missing values provided for some values. Check your inputs")
  }
  
  class(out) = 'simmr_input'
  
  return(out)
  
}
