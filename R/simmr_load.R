#' Load data into simmr and check for errors
#'
#' @param data A named list with the elements \code{mixtures}, \code{source_names}, \code{source_means}, and \code{source_sds}. Optionally \code{correction_means}, \code{correction_sds}, and \code{concentration_means} can be provided. \code{simmr_load} will provide a note if these are not used.
#' @param fixed_covars A matrix of covariates to be used as fixed effects. Should have same number of rows as \code{data$mixtures}.
#' @param random_covars A matrix of covariates to be used as random effects. Should have same number of rows as \code{data$mixtures}.
#' @param time_var A vector of time values (may be discrete or continuous) to use in a time series analysis. Should the same length as the number of rows as \code{data$mixtures}.
#'
#' @return A named list which is suitable for input into \link{\code{plot.simmr_input}}, \link{\code{simmr_mcmc}} and other fitting functions.
#' 
#' @export
#'
#' import tidyverse
#'
#' @examples
#' #' # A simple example with 10 observations, 2 tracers and 4 sources
#' 
#' data = list(
#' mixtures = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54,
#' -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89,
#' 11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10),
#' source_names=c('Source A','Source B','Source C','Source D'),
#' tracer_names = c('d13C','d15N'),
#' source_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, nrow=4),
#' source_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4),
#' correction_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol=2, nrow=4),
#' correction_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol=2, nrow=4),
#' concentration_means = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4))
#' 
#' # Load in with simmr_load
#' simmr_1 = data %>% simmr_load
#'
#' print(simmr_1)
#' 
simmr_load = function(data,
                      fixed_covars = NULL,
                      random_covars = NULL,
                      time_var = NULL,
                      print = TRUE,
                      iso_plot = TRUE) {
  
# New function to load in data for simmr so that it can be run using tidyverse type code
  
# Data needs to have at the very least the following elements
mix = data$mixtures
if(is.null(mix)) stop("Mixtures element of data is NULL")
s_names = data$source_names
if(is.null(s_names)) stop("source_names element of data is NULL")
tracer_names = data$tracer_names
if(is.null(tracer_names)) tracer_names = c('tracer_1', 'tracer_2')
s_means = data$source_means
if(is.null(s_means)) stop("source_means element of data is NULL")
s_sds = data$source_sds
if(is.null(s_sds)) stop("source_sds element of data is NULL")
c_means = data$correction_means
if(is.null(c_means)) cat("NOTE: correction_means element of data is NULL")
c_sds = data$correction_sds
if(is.null(c_sds)) cat("NOTE: correction_sds element of data is NULL")
conc_means = data$concentration_means
if(is.null(conc_means)) cat("NOTE: concentration_means element of data is NULL")

if (length(dim(mix)) != 2)
  stop("mixtures object must be a 2D object")

n_obs = dim(mix)[1]
n_tracers = dim(mix)[2]

# Add column names if they're not there
if (is.null(colnames(mix)))
  colnames(mix) = paste0('tracer', 1:n_tracers)

# source_names must be a character vector - the length of it is the number of sources
if (!is.vector(s_names))
  stop("source_names must be a vector")
s_names = as.character(s_names)
n_sources = length(s_names)

# source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
if (length(dim(s_means)) != 2)
  stop("source_means must have two dimensions")
if (length(dim(s_sds)) != 2)
  stop("source_sds must have two dimensions")
if (nrow(s_means) != n_sources)
  stop('Number of rows in source_means does not match length(source_names)')
if (ncol(s_means) != n_tracers)
  stop('Number of columns in source_means does not match ncol(mixtures)')
if (nrow(s_sds) != n_sources)
  stop('Number of rows in source_sds does not match length(source_names)')
if (ncol(s_sds) != n_tracers)
  stop('Number of columns in source_sds does not match ncol(mixtures)')

# Check that either neither or both of corrections_means and sds are given
if (is.null(c_means) &
    !is.null(c_sds))
  stop("Both correction_means and correction_sds must be supplied")
if (!is.null(c_means) &
    is.null(c_sds))
  stop("Both correction_means and correction_sds must be supplied")

# correction_means and correction_sds must be matrix and of same dimension as n_sources if they supplied
if (!is.null(c_means)) {
  if (length(dim(c_means)) != 2)
    stop("correction_means must have two dimensions")
  if (length(dim(c_sds)) != 2)
    stop("correction_sds must have two dimensions")
  if (nrow(c_sds) != n_sources)
    stop('Number of rows in correction_means does not match length(source_names)')
  if (ncol(c_means) != n_tracers)
    stop('Number of columns in correction_means does not match ncol(mixtures)')
  if (nrow(c_sds) != n_sources)
    stop('Number of rows in correction_sds does not match length(source_names)')
  if (ncol(c_sds) != n_tracers)
    stop('Number of columns in correction_sds does not match ncol(mixtures)')
} else {
  c_means = matrix(0, ncol = n_tracers, nrow = n_sources)
  c_sds = matrix(0, ncol = n_tracers, nrow = n_sources)
}
  
  # concentration_means must be a matrix where all elements are less than 1
  if (!is.null(conc_means)) {
    if (length(dim(conc_means)) != 2)
      stop("concentration_means must have two dimensions")
    if (nrow(conc_means) != n_sources)
      stop('Number of rows in concentration_means does not match length(source_names)')
    if (ncol(conc_means) != n_tracers)
      stop('Number of columns in concentration_means does not match ncol(mixtures)')
    if (any(conc_means > 1))
      stop("Invalid values in concentration_means; must all be less than 1")
    if (any(conc_means == 0))
      stop("Invalid values in concentration_means; must all be non-zero")
    
  } else {
    conc_means = matrix(1, ncol = n_tracers, nrow = n_sources)
  }
  
  # Sort out fixed factors
  # fixed_covars
  fixed_mat = NULL
  
  # Sort out random factors
  random_mat = NULL
  
  # Sort out time variable
  time_var = NULL

  # Prepare output and give class
  out = list(
    mixtures = mix,
    tracer_names = tracer_names,
    source_names = s_names,
    source_means = s_means,
    source_sds = s_sds,
    correction_means = c_means,
    correction_sds = c_sds,
    concentration_means = conc_means,
    fixed_mat = fixed_mat,
    random_mat = random_mat,
    time_var = time_var,
    n_obs = n_obs,
    n_tracers = n_tracers,
    n_sources = n_sources
  )
  
  # Look through to see whether there are any missing values in anything
  if (any(unlist(lapply(out, 'is.na')))) {
    warning("Missing values provided for some elements Check your inputs")
  }
  
  # Put into the correct class
  if(!is.null(fixed_mat) | !is.null(random_mat)) {
    class(out) = 'simmr_input_re'
  } else if(!is.null(time_var)) {
    class(out) = 'simmr_input_ts'
  } else {
    class(out) = 'simmr_input'  
  }
  
  if(print) print.simmr_input(out)
  if(iso_plot) plot.simmr_input(out)
  
  return(out)
  
}
