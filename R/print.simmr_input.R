#' Print simmr input object
#'
#' @param x An object of class \code{simmr_input}
#' @param ... Other arguments (not supported)
#'
#' @return A neat presentation of your simmr object.
#' @export
print.simmr_input <-
  function(x, ...) {
    message("This is a valid simmr input object with ")
    message(paste(x$n_obs, " observations, "))
    message(paste(x$n_tracers, "tracers, and "))
    message(paste(x$n_sources, "sources.\n"))
    if (x$n_groups > 1) cat(paste("There are", x$n_groups, "groups.\n"))
    message("The source names are: ")
    print(x$source_names, sep = ", ")
    message(".\n")
    message("The tracer names are: ")
    print(colnames(x$mixtures), sep = ", ")
    message("\n\n")
  }
