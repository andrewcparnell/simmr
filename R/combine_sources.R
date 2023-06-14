#' Combine the dietary proportions from two food sources after running simmr
#'
#' This function takes in an object of class \code{simmr_output} and combines
#' two of the food sources. It works for single and multiple group data.
#'
#' Often two sources either (1) lie in similar location on the iso-space plot,
#' or (2) are very similar in phylogenetic terms. In case (1) it is common to
#' experience high (negative) posterior correlations between the sources.
#' Combining them can reduce this correlation and improve precision of the
#' estimates. In case (2) we might wish to determine the joint amount eaten of
#' the two sources when combined. This function thus combines two sources after
#' a run of \code{\link{simmr_mcmc}} or \code{\link{simmr_ffvb}} (known as
#' a posteriori combination). The new object can then be called with
#' \code{\link{plot.simmr_input}} or \code{\link{plot.simmr_output}} to
#' produce iso-space plots of summaries of the output after combination.
#'
#' @param simmr_out An object of class \code{simmr_output} created from
#' \code{\link{simmr_mcmc}} or \code{\link{simmr_ffvb}}
#' @param to_combine The names of exactly two sources. These should match the
#' names given to \code{\link{simmr_load}}.
#' @param new_source_name A name to give to the new combined source.
#' @return A new \code{simmr_output} object
#' @author Andrew Parnell <andrew.parnell@@mu.ie>, Emma Govan
#' @seealso See \code{\link{simmr_mcmc}} and \code{\link{simmr_ffvb}} and
#' the associated vignette for examples.
#' @examples
#' \donttest{
#' # The data
#' data(geese_data)
#'
#' # Load into simmr
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
#' # Plot
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Print it
#' print(simmr_1_out)
#'
#' # Summary
#' summary(simmr_1_out)
#' summary(simmr_1_out, type = "diagnostics")
#' summary(simmr_1_out, type = "correlations")
#' summary(simmr_1_out, type = "statistics")
#' ans <- summary(simmr_1_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(simmr_1_out)
#' plot(simmr_1_out, type = "boxplot")
#' plot(simmr_1_out, type = "histogram")
#' plot(simmr_1_out, type = "density")
#' plot(simmr_1_out, type = "matrix")
#'
# simmr_out_combine <- combine_sources(simmr_1_out,
#   to_combine = c("U.lactuca", "Enteromorpha"),
#   new_source_name = "U.lac+Ent"
# )
#' plot(simmr_out_combine$input)
# plot(simmr_out_combine, type = "boxplot", title = "simmr output: combined sources")
# }
#' @export
combine_sources <- function(simmr_out,
                            to_combine = simmr_out$input$source_names[1:2],
                            new_source_name = "combined_source") {
  UseMethod("combine_sources")
}
#' @export
combine_sources.simmr_output <- function(simmr_out,
                                         to_combine = simmr_out$input$source_names[1:2],
                                         new_source_name = "combined_source") {
  if (inherits(simmr_out, "mcmc") == TRUE) {
    # Check that to_combine is in the list of sources
    assert_true(all(to_combine %in% simmr_out$input$source_names))

    # Check only two sources to be combined
    assert_character(to_combine,
      any.missing = FALSE,
      all.missing = FALSE,
      max.len = simmr_out$input$n_sources - 1,
      unique = TRUE
    )

    # If there are more than two sources call the function recursively
    if (length(to_combine) > 2) {
      curr_new_source_name <- paste0(to_combine[-1], collapse = "+")
      simmr_out <- combine_sources(simmr_out,
        to_combine = to_combine[-1],
        new_source_name = curr_new_source_name
      )
      to_combine <- c(curr_new_source_name, to_combine[1])
    }

    # Find which columns to combine by number
    to_combine_cols <- sort(match(to_combine, simmr_out$input$source_names))

    # Create a new object
    simmr_new_out <- simmr_out

    # 1 combine the chosen source means
    old_source_means <- simmr_out$input$source_means
    simmr_new_out$input$source_means <- old_source_means[-to_combine_cols[2], ,
      drop = FALSE
    ]
    simmr_new_out$input$source_means[to_combine_cols[1], ] <- matrix(apply(old_source_means[to_combine_cols, , drop = FALSE], 2, "mean"), nrow = 1)

    # 2 combine the source sds
    old_source_sds <- simmr_out$input$source_sds
    simmr_new_out$input$source_sds <- old_source_sds[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$source_sds[to_combine_cols[1], ] <- matrix(apply(old_source_sds[to_combine_cols, , drop = FALSE], 2, function(x) {
      sqrt(sum(x^2))
    }), nrow = 1)

    # 3 combine the correction means
    old_correction_means <- simmr_out$input$correction_means
    simmr_new_out$input$correction_means <- old_correction_means[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$correction_means[to_combine_cols[1], ] <- matrix(apply(old_correction_means[to_combine_cols, , drop = FALSE], 2, "mean"), nrow = 1)

    # 4 combine the correction sds
    old_correction_sds <- simmr_out$input$correction_sds
    simmr_new_out$input$correction_sds <- old_correction_sds[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$correction_sds[to_combine_cols[1], ] <- matrix(apply(old_correction_sds[to_combine_cols, , drop = FALSE], 2, function(x) {
      sqrt(sum(x^2))
    }), nrow = 1)

    # 5 combine the concentration means
    old_concentration_means <- simmr_out$input$concentration_means
    simmr_new_out$input$concentration_means <- old_concentration_means[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$concentration_means[to_combine_cols[1], ] <- matrix(apply(old_concentration_means[to_combine_cols, , drop = FALSE], 2, "mean"), nrow = 1)

    # 6 change the source names
    old_source_names <- simmr_out$input$source_names
    simmr_new_out$input$source_names <- old_source_names[-to_combine_cols[2]]
    simmr_new_out$input$source_names[to_combine_cols[1]] <- new_source_name

    # 7 Change n_sources
    simmr_new_out$input$n_sources <- simmr_new_out$input$n_sources - 1



    # 8 Sum across all the output values
    n_groups <- simmr_out$input$n_groups
    for (j in 1:n_groups) {
      simmr_new_out$output[[j]] <- simmr_out$output[[j]]
      # Change sims.list and sims.matrix
      # First sims.matrix
      sims.matrix <- simmr_out$output[[j]]$BUGSoutput$sims.matrix
      new_sims.matrix <- sims.matrix[, -(to_combine_cols[2] + 1)]
      new_sims.matrix[, to_combine_cols[1] + 1] <- sims.matrix[, to_combine_cols[1] + 1] + sims.matrix[, to_combine_cols[2] + 1]
      colnames(new_sims.matrix)[to_combine_cols[1] + 1] <- new_source_name
      simmr_new_out$output[[j]]$BUGSoutput$sims.matrix <- new_sims.matrix
      # Now sims.list
      sims.list <- simmr_out$output[[j]]$BUGSoutput$sims.list
      new_sims.list <- sims.list
      new_sims.list$p <- sims.list$p[, -to_combine_cols[2]]
      new_sims.list$p[, to_combine_cols[1]] <- sims.list$p[, to_combine_cols[2]] + sims.list$p[, to_combine_cols[1]]
      colnames(new_sims.list$p)[to_combine_cols[1]] <- new_source_name
      simmr_new_out$output[[j]]$BUGSoutput$sims.list <- new_sims.list
    }
  }
  if (inherits(simmr_out, "ffvb") == TRUE) {
    # Check that to_combine is in the list of sources
    assert_true(all(to_combine %in% simmr_out$input$source_names))

    # Check only two sources to be combined
    assert_character(to_combine,
      any.missing = FALSE,
      all.missing = FALSE,
      max.len = simmr_out$input$n_sources - 1,
      unique = TRUE
    )

    # If there are more than two sources call the function recursively
    if (length(to_combine) > 2) {
      curr_new_source_name <- paste0(to_combine[-1], collapse = "+")
      simmr_out <- combine_sources(simmr_out,
        to_combine = to_combine[-1],
        new_source_name = curr_new_source_name
      )
      to_combine <- c(curr_new_source_name, to_combine[1])
    }

    # Find which columns to combine by number
    to_combine_cols <- sort(match(to_combine, simmr_out$input$source_names))

    # Create a new object
    simmr_new_out <- simmr_out

    # 1 combine the chosen source means
    old_source_means <- simmr_out$input$source_means
    simmr_new_out$input$source_means <- old_source_means[-to_combine_cols[2], ,
      drop = FALSE
    ]
    simmr_new_out$input$source_means[to_combine_cols[1], ] <- matrix(apply(old_source_means[to_combine_cols, , drop = FALSE], 2, "mean"), nrow = 1)

    # 2 combine the source sds
    old_source_sds <- simmr_out$input$source_sds
    simmr_new_out$input$source_sds <- old_source_sds[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$source_sds[to_combine_cols[1], ] <- matrix(apply(old_source_sds[to_combine_cols, , drop = FALSE], 2, function(x) {
      sqrt(sum(x^2))
    }), nrow = 1)

    # 3 combine the correction means
    old_correction_means <- simmr_out$input$correction_means
    simmr_new_out$input$correction_means <- old_correction_means[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$correction_means[to_combine_cols[1], ] <- matrix(apply(old_correction_means[to_combine_cols, , drop = FALSE], 2, "mean"), nrow = 1)

    # 4 combine the correction sds
    old_correction_sds <- simmr_out$input$correction_sds
    simmr_new_out$input$correction_sds <- old_correction_sds[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$correction_sds[to_combine_cols[1], ] <- matrix(apply(old_correction_sds[to_combine_cols, , drop = FALSE], 2, function(x) {
      sqrt(sum(x^2))
    }), nrow = 1)

    # 5 combine the concentration means
    old_concentration_means <- simmr_out$input$concentration_means
    simmr_new_out$input$concentration_means <- old_concentration_means[-to_combine_cols[2], , drop = FALSE]
    simmr_new_out$input$concentration_means[to_combine_cols[1], ] <- matrix(apply(old_concentration_means[to_combine_cols, , drop = FALSE], 2, "mean"), nrow = 1)

    # 6 change the source names
    old_source_names <- simmr_out$input$source_names
    simmr_new_out$input$source_names <- old_source_names[-to_combine_cols[2]]
    simmr_new_out$input$source_names[to_combine_cols[1]] <- new_source_name

    # 7 Change n_sources
    simmr_new_out$input$n_sources <- simmr_new_out$input$n_sources - 1

    # 8 Sum across all the output values
    n_groups <- simmr_out$input$n_groups
    for (j in 1:n_groups) {
      simmr_new_out$output[[j]] <- simmr_out$output[[j]]
      # Change sims.list and sims.matrix
      # First sims.matrix
      sims.matrix <- simmr_out$output[[j]]$BUGSoutput$sims.matrix
      new_sims.matrix <- sims.matrix[, -(to_combine_cols[2])]
      new_sims.matrix[, to_combine_cols[1]] <- sims.matrix[, to_combine_cols[1]] + sims.matrix[, to_combine_cols[2]]
      colnames(new_sims.matrix)[to_combine_cols[1]] <- new_source_name
      simmr_new_out$output[[j]]$BUGSoutput$sims.matrix <- new_sims.matrix
      # Now sims.list
      sims.list <- simmr_out$output[[j]]$BUGSoutput$sims.list
      new_sims.list <- sims.list
      new_sims.list$p <- sims.list$p[, -to_combine_cols[2]]
      new_sims.list$p[, to_combine_cols[1]] <- sims.list$p[, to_combine_cols[2]] + sims.list$p[, to_combine_cols[1]]
      colnames(new_sims.list$p)[to_combine_cols[1]] <- new_source_name
      simmr_new_out$output[[j]]$BUGSoutput$sims.list <- new_sims.list
    }
  }


  a <- c("ffvb", "mcmc")
  if (inherits(simmr_out, "ffvb") == TRUE) {
    i <- 1
  } else if (inherits(simmr_out, "mcmc") == TRUE) {
    i <- 2
  }
  class(simmr_new_out) <- c("simmr_output", a[i])
  return(simmr_new_out)
}
