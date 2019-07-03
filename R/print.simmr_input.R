#' Print simmr input object
#'
#' @param x An object of class \code{simmr_input}
#' @param ... Other arguments (not supported)
#'
#' @return A neat presentation of your simmr object. 
#' @export
print.simmr_input <-
function(x,...) {
  cat('Valid simmr input object with:\n')
  cat(paste(x$n_obs,'observations, '))
  cat(paste(x$n_tracers,'tracers, and '))
  cat(paste(x$n_sources,'sources.\n'))
  if(!is.null(x$correction_means)) {
    cat('It contains correction means and sds.\n')
  } else {
    cat('It does not contain correction means or sds.\n')
  }
  if(!is.null(x$concentration_means)) {
    cat('It also contains concentration means.\n')  
  } else {
    cat('It does not contain concentration means.\n')
  }
  
  
  cat('The source names are: ')
  cat(x$source_names,sep=', ')
  cat('.\n')
  cat('The tracer names are: ')
  cat(x$tracer_names, sep = ', ')
  cat('.\n')
}