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
#' @import checkmate
#'
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
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' print(simmr_1)
#' @export simmr_load
simmr_load <- function(mixtures,
                       source_names,
                       source_means,
                       source_sds,
                       correction_means = NULL,
                       correction_sds = NULL,
                       concentration_means = NULL,
                       group = NULL) {
  # Function to load in data for simmr and check whether it's appropriate for running through simmr_mcmc

  # Go through each object and check that it matches the requirements

  # Write a function that generically tests for any 2D numeric data shape such as matrix, data frame or tibble
  assert_2D_numeric <- function(x,
                                nrows = NULL,
                                ncols = NULL,
                                null.ok = FALSE) {
    assert(
      test_data_frame(x,
        types = c("double", "numeric"),
        nrows = nrows,
        ncols = ncols,
        null.ok = null.ok
      ),
      test_matrix(x,
        mode = "numeric",
        nrows = nrows,
        ncols = ncols,
        null.ok = null.ok
      ),
      test_tibble(x,
        types = c("double", "numeric"),
        nrows = nrows,
        ncols = ncols,
        null.ok = null.ok
      )
    )
  }

  # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
  # assert_matrix(mixtures)
  assert_2D_numeric(mixtures)
  n_obs <- nrow(mixtures)
  n_tracers <- ncol(mixtures)

  # Add column names if they're not there
  if (is.null(colnames(mixtures))) {
    colnames(mixtures) <- paste0("tracer", 1:n_tracers)
  }

  # source_names must be a character vector - the length of it is the number of sources
  assert_character(source_names)
  n_sources <- length(source_names)

  # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
  assert_2D_numeric(source_means,
    nrows = n_sources,
    ncols = n_tracers
  )
  # assert_matrix(source_means, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(source_sds,
    nrows = n_sources,
    ncols = n_tracers
  )
  # assert_matrix(source_sds, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(correction_means,
    nrows = n_sources,
    ncols = n_tracers,
    null.ok = ifelse(is.null(correction_sds),
      TRUE, FALSE
    )
  )
  # assert_matrix(correction_means,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_sds),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(correction_sds,
    nrows = n_sources,
    ncols = n_tracers,
    null.ok = ifelse(is.null(correction_sds),
      TRUE, FALSE
    )
  )
  # assert_matrix(correction_sds,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_means),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(concentration_means,
    nrows = n_sources,
    ncols = n_tracers,
    null.ok = TRUE
  )
  # assert_matrix(concentration_means,
  #   nrows = n_sources,
  #   ncols = n_tracers, null.ok = TRUE
  # )

  # Fill in correction means
  if (is.null(correction_means)) {
    correction_means <- matrix(0, ncol = n_tracers, nrow = n_sources)
    correction_sds <- matrix(0, ncol = n_tracers, nrow = n_sources)
  }

  # concentration_means must be a matrix where all elements are less than 1
  if (is.null(concentration_means)) {
    concentration_means <- matrix(1, ncol = n_tracers, nrow = n_sources)
  } else {
    assert_true(all(concentration_means < 1) & all(concentration_means > 0))
  }

  # Check the groups are the right length and structure if given
  if (!is.null(group)) {
    # Allow for the group to be a factor variable
    group <- as.factor(group)
    group_int <- as.integer(group)
    assert_vector(group, len = n_obs)
  } else {
    group <- factor(rep(1, n_obs))
    group_int <- as.integer(group)
  }
  n_groups <- length(unique(group))

  # Prepare output and give class
  out <- list(
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
  if (any(unlist(lapply(out, "is.na")))) {
    warning("Missing values provided for some values. Check your inputs")
  }

  class(out) <- "simmr_input"

  return(out)
}
