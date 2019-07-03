#' Print simmr input object
#'
#' @param x An object of class \code{simmr_input}
#' @param ... Other arguments (not supported)
#'
#' @return A neat presentation of your simmr object. 
#' @export
print.simmr_input <-
function(x,...) {
  cat('This is a valid simmr input object with ')
  cat(paste(x$n_obs,'observations, '))
  cat(paste(x$n_tracers,'tracers, and '))
  cat(paste(x$n_sources,'sources.\n'))
  if(x$n_groups>1) cat(paste('There are',x$n_groups,'groups.\n'))
  cat('The source names are: ')
  cat(x$source_names,sep=', ')
  cat('.\n')
  cat('The tracer names are: ')
  cat(colnames(x$mixtures), sep = ', ')
  cat('.\n\n')
}
