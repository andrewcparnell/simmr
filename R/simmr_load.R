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
#' @importFrom dplyr '%>%'
#'
#' @examples
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # See simmr_mcmc and vignettes for full example run
#' data(simmr_data_1)
#'
#' # Load into simmr
#' simmr_1 = simmr_data_1 %>% 
#'   simmr_load
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

  # First part is error checking
  
  # Data needs to have at the very least the following elements
  if(is.null(data$mixtures))
    stop("Mixtures element of data is missing")
  mix = as.matrix(data$mixtures)  
  n_obs = dim(mix)[1]
  n_tracers = dim(mix)[2]
  

  if(is.null(data$source_names))
    stop("source_names element of data is missing")
  s_names = as.character(unlist(data$source_names))
  n_sources = length(s_names)
  
  if(is.null(data$tracer_names)) {
    cat("NOTE: tracer_names element of data is missing Using tracer_X instead.\n")
    tracer_names = paste0('tracer_', 1:nrow(mix))
  } else {
    tracer_names = as.character(unlist(data$tracer_names))
  }
  
  if(is.null(data$source_means))
    stop("source_means element of data is missing")
  s_means = as.matrix(data$source_means)
  
  if(is.null(data$source_sds))
    stop("source_sds element of data is missing")
  s_sds = as.matrix(data$source_sds)
  
  if(is.null(data$correction_means)) {
    cat("NOTE: correction_means element of data is missing.\n")
    c_means = matrix(0, ncol = n_tracers, nrow = n_sources)
  } else {
    c_means = as.matrix(data$correction_means)
  }

  if(is.null(data$correction_sds)) {
    cat("NOTE: correction_sds element of data is missing.\n")
    c_sds = matrix(0, ncol = n_tracers, nrow = n_sources)
  } else {
    c_sds = as.matrix(data$correction_sds)
  }
  
  if(is.null(data$concetration_means)) {
    cat("NOTE: concetration_means element of data is missing.\n")
    conc_means = matrix(1, ncol = n_tracers, nrow = n_sources)
  } else {
    conc_means = as.matrix(data$concetration_means)
  }

if(length(dim(mix)) != 2)
  stop("mixtures object must be a 2D object")

# Add column names if they're not there
if(is.null(colnames(mix)))
  colnames(mix) = paste0('tracer', 1:n_tracers)

# source_names must be a character vector - the length of it is the number of sources
if (!is.vector(s_names))
  stop("source_names must be a vector")

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
  
  if(print) print(out)
  if(iso_plot) plot(out)
  
  return(out)
  
}
