#' Compare dietary proportions for a single source across different groups
#'
#' This function takes in an object of class \code{simmr_output} and creates
#' probabilistic comparisons for a given source and a set of at least two
#' groups.
#'
#' When two groups are specified, the function produces a direct calculation of
#' the probability that one group is bigger than the other. When more than two
#' groups are given, the function produces a set of most likely probabilistic
#' orderings for each combination of groups. The function produces boxplots by
#' default and also allows for the storage of the output for further analysis if
#' required.
#'
#' @param simmr_out An object of class \code{simmr_output} created from
#' \code{\link{simmr_mcmc}}.
#' @param source_name The name of a source. This should match the names exactly
#' given to \code{\link{simmr_load}}.
#' @param groups The integer values of the group numbers to be compared. At
#' least two groups must be specified.
#' @param plot A logical value specifying whether plots should be produced or
#' not.
#'
#' @import ggplot2
#'
#' @return If there are two groups, a vector containing the differences between
#' the two groups proportions for that source. If there are multiple groups, a
#' list containing the following fields: \item{Ordering }{The different
#' possible orderings of the dietary proportions across groups} \item{out_all
#' }{The dietary proportions for this source across the groups specified as
#' columns in a matrix}
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} for complete examples.
#' @examples
#' \dontrun{
#' data(geese_data)
#' simmr_in <- with(
#'   geese_data,
#'   simmr_load(
#'     mixtures = mixtures,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means,
#'     group = groups
#'   )
#' )
#'
#' # Print
#' simmr_in
#'
#' # Plot
#' plot(simmr_in,
#'   group = 1:8, xlab = expression(paste(delta^13, "C (\\u2030)", sep = "")),
#'   ylab = expression(paste(delta^15, "N (\\u2030)", sep = "")),
#'   title = "Isospace plot of Inger et al Geese data"
#' )
#'
#' # Run MCMC for each group
#' simmr_out <- simmr_mcmc(simmr_in)
#'
#' # Print output
#' simmr_out
#'
#' # Summarise output
#' summary(simmr_out, type = "quantiles", group = 1)
#' summary(simmr_out, type = "quantiles", group = c(1, 3))
#' summary(simmr_out, type = c("quantiles", "statistics"), group = c(1, 3))
#'
#' # Plot - only a single group allowed
#' plot(simmr_out, type = "boxplot", group = 2, title = "simmr output group 2")
#' plot(simmr_out, type = c("density", "matrix"), grp = 6, title = "simmr output group 6")
#'
#' # Compare groups
#' compare_groups(simmr_out, source = "Zostera", groups = 1:2)
#' compare_groups(simmr_out, source = "Zostera", groups = 1:3)
#' compare_groups(simmr_out, source = "U.lactuca", groups = c(4:5, 7, 2))
#' }
#'
#' @export
compare_groups <- function(simmr_out,
                           source_name = simmr_out$input$source_names[1],
                           groups = 1:2,
                           plot = TRUE) {
  UseMethod("compare_groups")
}
#' @export
compare_groups.simmr_output <- function(simmr_out,
                                        source_name = simmr_out$input$source_names[1],
                                        groups = 1:2,
                                        plot = TRUE) {

  # Function to compare between groups both via textual output and with boxplots
  # Things to supply are:
  # If two groups are given:
  #   - provide the probability of one group being bigger than the other
  #   - give the probability distribution of the difference
  #   - optional boxplot of two
  # If more than two groups are given:
  #   - provide the top most likely orderings of the groups
  # An optional boxplot of the groups

  # Throw an error if only one group is specified
  assert_true(simmr_out$input$n_groups > 1)
  assert_numeric(groups, min.len = 2)
  assert_true(all(source_name %in% simmr_out$input$source_names))
  assert_character(source_name,
    any.missing = FALSE,
    len = 1
  )
  assert_true(all(source_name %in% simmr_out$input$source_names))

  # Get group names
  group_names <- levels(simmr_out$input$group)[groups]

  # Start with two groups version
  if (length(groups) == 2) {
    # Get the output for this particular source on these two groups
    post_mat1 <- simmr_out$output[[groups[1]]]$BUGSoutput$sims.matrix
    post_mat2 <- simmr_out$output[[groups[2]]]$BUGSoutput$sims.matrix
    match_name <- match(source_name, colnames(post_mat1))
    out_all_grp_1 <- post_mat1[, match_name]
    out_all_grp_2 <- post_mat2[, match_name]
    # Produce the difference between the two
    out_diff <- out_all_grp_1 - out_all_grp_2

    cat(paste("Prob ( proportion of", source_name, "in group", group_names[1], "> proportion of", source_name, "in group", group_names[2], ") =", round(mean(out_diff > 0), 3)))

    if (plot) {
      # Stupid fix for packaging ggplot things
      Group <- Proportion <- NULL
      df <- data.frame(Proportion = c(out_all_grp_1, out_all_grp_2), Group = c(rep(paste(group_names[1]), length(out_all_grp_1)), rep(paste(group_names[2]), length(out_all_grp_2))))
      p <- ggplot(df, aes(x = Group, y = Proportion, fill = Group)) +
        geom_boxplot(alpha = 0.5, outlier.size = 0) +
        theme_bw() +
        theme(legend.position = "none") +
        ggtitle(paste("Comparison of dietary proportions for groups", group_names[1], "and", group_names[2], "for source", source_name))
      print(p)
    }
  }

  # Now for more groups
  if (length(groups) > 2) {
    # Get the output for all the groups
    post_mat <- simmr_out$output[[groups[1]]]$BUGSoutput$sims.matrix
    match_name <- match(source_name, colnames(post_mat))
    len <- length(post_mat[, match_name])
    out_all <- matrix(NA, nrow = len, ncol = length(groups))
    for (j in 1:length(groups)) out_all[, j] <- simmr_out$output[[groups[j]]]$BUGSoutput$sims.matrix[, match_name]
    colnames(out_all) <- paste(group_names)

    # Now find the ordering of each one
    ordering_num <- t(apply(out_all, 1, order, decreasing = TRUE))
    Ordering <- rep(NA, length = nrow(ordering_num))
    for (i in 1:length(Ordering)) Ordering[i] <- paste0(group_names[ordering_num[i, ]], collapse = " > ")
    cat("Most popular orderings are as follows:\n")
    tab <- t(t(sort(table(Ordering, dnn = NULL), decreasing = TRUE)))
    colnames(tab) <- "Probability"
    # Do not print all of it if too long
    if (nrow(tab) > 30) {
      print(round(tab[1:30, ] / length(Ordering), 4))
    } else {
      print(round(tab / length(Ordering), 4))
    }

    if (plot) {
      # Stupid fix for packaging ggplot things
      Group <- Proportion <- NULL
      df <- reshape2::melt(out_all)[, 2:3]
      colnames(df) <- c("Group", "Proportion")
      p <- ggplot(df, aes(x = factor(Group), y = Proportion, fill = factor(Group))) +
        scale_fill_viridis(discrete = TRUE) +
        geom_boxplot(alpha = 0.5, outlier.size = 0) +
        xlab("Group") +
        theme_bw() +
        theme(legend.position = "none") +
        ggtitle(paste("Comparison of dietary proportions for source", source_name))
      print(p)
    }
  }

  # Return output
  if (length(groups) == 2) {
    if (plot) {
      invisible(list(out_diff = out_diff, plot = p))
    } else {
      invisible(list(out_diff = out_diff))
    }
  } else {
    if (plot) {
      invisible(list(Ordering = Ordering, out_all = out_all, plot = p))
    } else {
      invisible(list(Ordering = Ordering, out_all = out_all))
    }
  }
}
