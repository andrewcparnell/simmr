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
#'
#'   @seealso \code{\link{simmr_mcmc}} for an example where it is used
#'
#' }
"simmr_data_1"

#' A 3-isotope fake stable isotope mixing data set
#'
#' A fake data set with 30 observations on 3 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A three column matrix containing delta 13C, delta 15N, and delta 34S values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N, d34S)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'   ...
#'
#'   @seealso \code{\link{simmr_mcmc}} for an example where it is used
#'
#' }
"simmr_data_2"

#' Geese stable isotope mixing data set
#'
#' A real Geese data set with 251 observations on 2 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means. Taken from Inger et al (2016). See link for paper
#'
#' @source \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2656.2006.01142.x}
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N, d34S)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'   ...
#'
#'   @seealso \code{\link{simmr_mcmc}} for an example where it is used
#' }
"geese_data"

#' An artificial data set used to indicate effect of priors
#'
#' A fake box data set identified by Fry (2014) as a failing of SIMMs
#' See the link for more interpretation of these data and the output
#'
#' @source \url{https://www.int-res.com/abstracts/meps/v490/p285-289/}
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'   ...
#'
#'   @seealso \code{\link{simmr_mcmc}} for an example where it is used
#'
#'
#' }
"square_data"

#' A smaller version of the Geese stable isotope mixing data set
#'
#' A real Geese data set with 9 observations on 2 isotopes, with 4 sources, and with corrections/trophic enrichment factors (TEFs or TDFs), and concentration dependence means. Taken from Inger et al (2016). See link for paper
#'
#' @source \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2656.2006.01142.x}
#'
#' @format A list with the following elements
#' \describe{
#'   \item{mixtures}{A two column matrix containing delta 13C and delta 15N values respectively}
#'   \item{source_names}{A character vector of the food source names}
#'   \item{tracer_names}{A character vector of the tracer names (d13C, d15N, d34S)}
#'   \item{source_means}{A matrix of source mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{source_sds}{A matrix of source sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_means}{A matrix of TEFs mean values for the tracers in the same order as \code{mixtures} above}
#'   \item{correction_sds}{A matrix of TEFs sd values for the tracers in the same order as \code{mixtures} above}
#'   \item{concentration_means}{A matrix of concentration dependence mean values for the tracers in the same order as \code{mixtures} above}
#'   ...
#'
#'   @seealso \code{\link{simmr_mcmc}} for an example where it is used
#' }
"geese_data_day1"
