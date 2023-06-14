#' Plot different features of an object created from \code{\link{simmr_mcmc}}
#' or \code{\link{simmr_ffvb}}.
#'
#' This function allows for 4 different types of plots of the simmr output
#' created from \code{\link{simmr_mcmc}} or \code{\link{simmr_ffvb}}. The
#' types are: histogram, kernel density plot, matrix plot (most useful) and
#' boxplot. There are some minor customisation options.
#'
#' The matrix plot should form a necessary part of any SIMM analysis since it
#' allows the user to judge which sources are identifiable by the model.
#' Further detail about these plots is provided in the vignette.
#' Some code from
#' https://stackoverflow.com/questions/14711550/is-there-a-way-to-change-the-color-palette-for-ggallyggpairs-using-ggplot
#' accessed March 2023
#'
#' @param x An object of class \code{simmr_output} created via
#' \code{\link{simmr_mcmc}} or \code{\link{simmr_ffvb}}.
#' @param type The type of plot required. Can be one or more of 'histogram',
#' 'density', 'matrix', or 'boxplot'
#' @param group Which group(s) to plot.
#' @param binwidth The width of the bins for the histogram. Defaults to 0.05
#' @param alpha The degree of transparency of the plots. Not relevant for
#' matrix plots
#' @param title The title of the plot.
#' @param ggargs Extra arguments to be included in the ggplot (e.g. axis limits)
#' @param ...  Currently not used
#'
#' @import ggplot2
#' @import graphics
#' @import viridis
#' @importFrom reshape2 "melt"
#' @importFrom stats "cor"
#'
#' @author Andrew Parnell <andrew.parnell@@mu.ie>, Emma Govan
#' @seealso See \code{\link{simmr_mcmc}} and \code{\link{simmr_ffvb}} for
#' creating objects suitable for this function, and many more examples. See
#' also \code{\link{simmr_load}} for creating simmr objects,
#' \code{\link{plot.simmr_input}} for creating isospace plots,
#' \code{\link{summary.simmr_output}} for summarising output.
#'
#' @examples
#' \donttest{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
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
#' # Plot
#' plot(simmr_1)
#'
#'
#' # MCMC run
#' simmr_1_out <- simmr_mcmc(simmr_1)
#'
#' # Plot
#' plot(simmr_1_out) # Creates all 4 plots
#' plot(simmr_1_out, type = "boxplot")
#' plot(simmr_1_out, type = "histogram")
#' plot(simmr_1_out, type = "density")
#' plot(simmr_1_out, type = "matrix")
#' }
#' @export
plot.simmr_output <-
  function(x,
           type = c(
             "isospace",
             "histogram",
             "density",
             "matrix",
             "boxplot"
           ),
           group = 1,
           binwidth = 0.05,
           alpha = 0.5,
           title = if (length(group) == 1) {
             "simmr output plot"
           } else {
             paste("simmr output plot: group", group)
           },
           ggargs = NULL,
           ...) {
    if (inherits(x, "simmr_output") == TRUE) {
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)

      # Iso-space plot is special as all groups go on one plot
      # Add in extra dots here as they can be sent to this plot function
      if ("isospace" %in% type) graphics::plot(x$input, group = group, title = title, ...)

      # Get group names
      group_names <- levels(x$input$group)[group]

      for (i in 1:length(group)) {
        # Prep data
        out_all <- x$output[[group[i]]]$BUGSoutput$sims.list$p
        colnames(out_all) <- x$input$source_names
        df <- reshape2::melt(out_all)
        colnames(df) <- c("Num", "Source", "Proportion")
        if ("histogram" %in% type) {
          g <- ggplot(df, aes(
            x = Proportion,
            fill = Source
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(aes(y = after_stat(density)), binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title[i]) +
            facet_wrap("~ Source") +
            theme(legend.position = "none") +
            ggargs
          print(g)
        }

        if ("density" %in% type) {
          g <- ggplot(df, aes(
            x = Proportion,
            fill = Source
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_density(aes(y = after_stat(density)), alpha = alpha, linetype = 0) +
            theme_bw() +
            theme(legend.position = "none") +
            ggtitle(title[i]) +
            ylab("Density") +
            facet_wrap("~ Source") +
            ggargs
          print(g)
        }

        if ("boxplot" %in% type) {
          g <- ggplot(df, aes(
            y = Proportion, x = Source,
            fill = Source, alpha = alpha
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_boxplot(alpha = alpha, notch = TRUE, outlier.size = 0) +
            theme_bw() +
            ggtitle(title[i]) +
            theme(legend.position = "none") +
            coord_flip() +
            ggargs
          print(g)
        }

        # if ('convergence'%in%type) {
        #   coda::gelman.plot(x$output[[group[i]]],transform=TRUE)
        # }

        if ("matrix" %in% type) {
          modified_bar <- function(data, mapping, ...) {
            GGally::ggally_barDiag(data, mapping, ..., fill = viridis(1), binwidth = 0.025) + coord_cartesian(xlim = c(0, 1)) + theme_bw()
          }
          modified_density <- function(data, mapping, ...) {
            ggplot(data = data, mapping = mapping, ...) +
              stat_density_2d(
                geom = "polygon", contour = TRUE,
                aes(fill = after_stat(..level..)),
                bins = 5,
              ) +
              scale_fill_viridis(alpha = 0.8) + 
              #scale_fill_distiller(palette = "Blues", direction = 1) +
              theme_bw() +
              scale_x_continuous(limits = c(0, 1)) +
              scale_y_continuous(limits = c(0, 1))
          }

          g <- GGally::ggpairs(data.frame(out_all),
            upper = list(continuous = GGally::wrap(modified_density)),
            diag = list(continuous = GGally::wrap(modified_bar)),
            lower = list(continuous = GGally::wrap("cor", stars = FALSE))
          )



          print(g)
        }
      }
      if (exists("g")) invisible(g)
    }
  }
