#' A simple fake stable isotope mixing data set
#'
#' A simple fake data set with 10 observations on 2 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C and d15N)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'   ...
#' }
"simmr_data_1"


