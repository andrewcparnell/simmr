#' Plot different features of an object created from \code{\link{simmr_mcmc}}.
#'
#' This function allows for 4 different types of plots of the simmr output
#' created from \code{\link{simmr_mcmc}}. The types are: histogram, kernel
#' density plot, matrix plot (most useful) and boxplot. There are some minor
#' customisation options.
#'
#' The matrix plot should form a necessary part of any SIMM analysis since it
#' allows the user to judge which sources are identifiable by the model.
#' Further detail about these plots is provided in the vignette.
#'
#' @param x An object of class \code{simmr_output} created via
#' \code{\link{simmr_mcmc}}
#' @param type The type of plot required. Can be one or more of 'histogram',
#' 'density', 'matrix', or 'boxplot'
<<<<<<< HEAD
#' @param group Which group(s) to plot.
=======
>>>>>>> 46d3aea0df59ff3ce395f948d5c25dd6198c8c82
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
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} for creating objects suitable for this
#' function, and many more examples. See also \code{\link{simmr_load}} for
#' creating simmr objects, \code{\link{plot.simmr_input}} for creating isospace
#' plots, \code{\link{summary.simmr_output}} for summarising output.
#' @examples
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # See simmr_mcmc and vignettes for full example run
#' data(simmr_data_1)
#' 
#' # Load into simmr and run
#' simmr_1 = simmr_data_1 %>% 
#'   simmr_load %>% 
#'   simmr_mcmc
#' 
#' # Plot output
#' simmr_1 %>% plot(type='boxplot')
#' simmr_1 %>% plot(type='histogram')
#' simmr_1 %>% plot(type='density')
#' simmr_1 %>% plot(type='matrix') # Often the most useful
#' @export
plot.simmr_output <-
function(x,
         type = c('isospace',
                  'histogram',
                  'density',
                  'matrix',
                  'boxplot'),
         binwidth = 0.05,
         alpha = 0.5,
         title = 'simmr output plot',
         ggargs = NULL,
          ...) {

  # Get the specified type
  type=match.arg(type,several.ok=TRUE)

  # Iso-space plot is special as all groups go on one plot
  # Add in extra dots here as they can be sent to this plot function
  if('isospace' %in% type) graphics::plot(x$input,title=title,...)

<<<<<<< HEAD
  # Get group names
  group_names = levels(x$input$group)[group]
  
  for(i in 1:length(group)) {

    # Prep data
    out_all = x$output[[group[i]]]$BUGSoutput$sims.list$p
    colnames(out_all) = x$input$source_names
    df = reshape2::melt(out_all)
    colnames(df) = c('Num','Source','Proportion')
    if ('histogram'%in%type) {
      g=ggplot(df,aes_string(x="Proportion",y="..density..",
                             fill="Source")) +
        scale_fill_viridis(discrete=TRUE) +
        geom_histogram(binwidth=binwidth,alpha=alpha) +
        theme_bw() +
        ggtitle(title[i]) +
        facet_wrap("~ Source") +
        theme(legend.position='none') +
        ggargs
      print(g)
    }
=======
  # Prep data
  out_all = x$output$BUGSoutput$sims.list$p
  colnames(out_all) = x$input$source_names
  df = reshape2::melt(out_all)
  colnames(df) = c('Num','Source','Proportion')
  if ('histogram'%in%type) {
    g=ggplot(df,aes_string(x="Proportion",y="..density..",
                           fill="Source")) +
      scale_fill_viridis(discrete=TRUE) +
      geom_histogram(binwidth=binwidth,alpha=alpha) +
      theme_bw() +
      ggtitle(title) +
      facet_wrap("~ Source") +
      theme(legend.position='none') +
      ggargs
    print(g)
  }

  if ('density'%in%type) {
    g=ggplot(df,aes_string(x="Proportion",y="..density..",
                           fill="Source")) +
      scale_fill_viridis(discrete=TRUE) +
      geom_density(alpha=alpha,linetype=0) +
      theme_bw() +
      theme(legend.position='none') +
      ggtitle(title)  +
      ylab("Density") +
      facet_wrap("~ Source") +
      ggargs
    print(g)
  }
>>>>>>> 46d3aea0df59ff3ce395f948d5c25dd6198c8c82

  if ('boxplot'%in%type) {
    g=ggplot(df,aes_string(y="Proportion",x="Source",
                           fill="Source",alpha="alpha")) +
      scale_fill_viridis(discrete=TRUE) +
      geom_boxplot(alpha=alpha,notch=TRUE,outlier.size=0) +
      theme_bw() +
      ggtitle(title) +
      theme(legend.position='none') +
      coord_flip() +
      ggargs
    print(g)
  }

  if ('matrix'%in%type) {
    # These taken from the help(pairs) file
    panel.hist <- function(x, ...) {
      usr <- graphics::par("usr"); on.exit(graphics::par(usr))
      graphics::par(usr = c(usr[1:2], 0, 1.5) )
      h <- graphics::hist(x, plot = FALSE)
      breaks <- h$breaks; nB <- length(breaks)
      y <- h$counts; y <- y/max(y)
      graphics::rect(breaks[-nB], 0, breaks[-1], y, col = "lightblue", ...)
    }
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
      usr <- graphics::par("usr"); on.exit(graphics::par(usr))
      graphics::par(usr = c(0, 1, 0, 1))
      r <- stats::cor(x, y)
      txt <- format(c(r, 0.123456789), digits = digits)[1]
      txt <- paste0(prefix, txt)
      if(missing(cex.cor)) cex.cor <- 0.8/graphics::strwidth(txt)
      graphics::text(0.5, 0.5, txt, cex = cex.cor * abs(r))
    }
    panel.contour <- function(x, y, ...)
    {
      usr <- graphics::par("usr"); on.exit(graphics::par(usr))
      graphics::par(usr = c(usr[1:2], 0, 1.5) )
      kd <- MASS::kde2d(x,y)
      kdmax <- max(kd$z)
      graphics::contour(kd,add=TRUE,drawlabels=FALSE,levels=c(kdmax*0.1,kdmax*0.25,kdmax*0.5,kdmax*0.75,kdmax*0.9))
    }
    graphics::pairs(out_all,xlim=c(0,1),ylim=c(0,1),main=title,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.contour)
  }

}
