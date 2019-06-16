#' Plot the \code{simmr_input} data created from \code{simmr_load}
#' 
#' This function creates iso-space (AKA tracer-space or delta-space) plots.
#' They are vital in determining whether the data are suitable for running in a
#' SIMM.
#' 
#' It is desirable to have the vast majority of the mixture observations to be
#' inside the convex hull defined by the food sources. When there are more than
#' two tracers (as in one of the examples below) it is recommended to plot all
#' the different pairs of the food sources. See the vignette for further
#' details of richer plots.
#' 
#' @param x An object created via the function \code{\link{simmr_load}}
#' @param tracers The choice of tracers to plot. If there are more than two
#' tracers, it is recommended to plot every pair of tracers to determine
#' whether the mixtures lie in the mixing polygon defined by the sources
#' @param title A title for the graph
#' @param xlab The x-axis label. By default this is assumed to be delta-13C but
#' can be made richer if required. See examples below.
#' @param ylab The y-axis label. By default this is assumed to be delta-15N in
#' per mil but can be changed as with the x-axis label
#' @param sigmas The number of standard deviations to plot on the source
#' values. Defaults to 1.
#' @param mix_name A optional string containing the name of the mixture
#' objects, e.g. Geese.
#' @param colour If TRUE (default) creates a plot. If not, puts the plot in
#' black and white
#' @param ggargs Extra arguments to be included in the ggplot (e.g. axis limits)
#' @param ...  Not used
#' 
#' @import ggplot2
#' @import viridis
#' 
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{plot.simmr_output}} for plotting the output of a
#' simmr run. See \code{\link{simmr_mcmc}} for running a simmr object once the
#' iso-space is deemed acceptable.
#' @examples
#' 
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # See simmr_mcmc and vignettes for full example run
#' data(simmr_data_1)
#'
#' # Load into simmr - automatically produces plot
#' simmr_1 = simmr_data_1 %>% 
#'   simmr_load %>%
#'   simmr_mcmc
#'   
#' # Plot it (again) from the output
#' simmr_1 %>% plot(type = 'isospace')
#'  
#' # Or plot it from the input
#' simmr_1$input %>% plot
#' 
#' ### A more complicated example with 30 obs, 3 tracers and 4 sources
#' data(simmr_data_2)
#' 
#' # Load in and plot by default the first two isotopes
#' simmr_2 = simmr_data_2 %>%
#'             simmr_load()
#' 
#' # To get at other isotopes call plot 3 times
#' simmr_2 %>% plot(type = 'isospace', tracers = c(1, 2))
#' # Now plot d15N vs d34S
#' simmr_2 %>% plot(type = 'isospace', tracers = c(2, 3))
#' # and finally d13C vs d34S
#' simmr_2 %>% plot(type = 'isospace', tracers = c(1, 3))
#' # See vignette('simmr') for fancier x-axis labels
#' 
#' @export
plot.simmr_input <-
function(x,
         tracers = c(1,2),
         title = 'Tracers plot',
         xlab = x$tracer_names[tracers[1]],
         ylab = x$tracer_names[tracers[2]], 
         sigmas = 1, 
         mix_name = 'Mixtures', 
         ggargs = NULL,
         colour = TRUE, 
         ...) {

# First get the mean corrected sources and the sd corrected sources
source_means_c = x$source_means + x$correction_means
source_sds_c = sqrt(x$source_sds^2 + x$correction_sds^2)

# Set up data frame for ggplot - have to do it this stupid way because of cran
x2=c(source_means_c[,tracers[1]],x$mixtures[,tracers[1]])
x_lower=c(source_means_c[,tracers[1]]-sigmas*source_sds_c[,tracers[1]],x$mixtures[,tracers[1]])
x_upper=c(source_means_c[,tracers[1]]+sigmas*source_sds_c[,tracers[1]],x$mixtures[,tracers[1]])

if(ncol(x$mixtures)>1) {
  y=c(source_means_c[,tracers[2]],x$mixtures[,tracers[2]])
  y_lower=c(source_means_c[,tracers[2]]-sigmas*source_sds_c[,tracers[2]],x$mixtures[,tracers[2]])
  y_upper=c(source_means_c[,tracers[2]]+sigmas*source_sds_c[,tracers[2]],x$mixtures[,tracers[2]])
}

Source=factor(c(x$source_names,
                rep(mix_name,nrow(x$mixtures))),
              levels=c(mix_name,x$source_names))
size=c(rep(0.5,x$n_sources),rep(0.5,nrow(x$mixtures)))

if(ncol(x$mixtures) == 1) {
  df=data.frame(x=x2,x_lower,x_upper,Source,size, y = Source)
} else {
  df=data.frame(x=x2,y=y,x_lower,y_lower,x_upper,y_upper,Source,size)
}

# Plot for bivariate mixtures
if(ncol(x$mixtures) > 1) {
  if(colour) {
    g=ggplot(data=df, aes(x = x,y = y,colour=Source)) + 
      scale_color_viridis(discrete=TRUE) + 
      theme_bw() +
      labs(x=xlab,y=ylab,title=title) +
      geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
      #geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
      geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,shape=Source)) +
      scale_shape_manual(values=1:nlevels(df$Source)) +
      theme(legend.title=element_blank(),legend.key = element_blank()) +
      guides(color=guide_legend(override.aes=list(linetype=c(rep(0,1),rep(1,x$n_sources))))) +
      ggargs
  } else {
    g=ggplot(data=df, aes(x = x,y = y,colour=Source)) + 
      theme_bw() +
      labs(x=xlab,y=ylab,title=title) +
      geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
      #geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
      geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,shape=Source)) +
      scale_shape_manual(values=1:nlevels(df$Source)) +
      theme(legend.title=element_blank(),legend.key = element_blank()) +
      guides(color=guide_legend(override.aes=list(linetype=c(rep(0,1),rep(1,x$n_sources))))) +
      scale_colour_grey() +
      ggargs
  }
}

# Plot for univariate mixtures
if(ncol(x$mixtures) == 1) {
  if(colour) {
    g = ggplot(data=df, aes(x = x, y = y, colour = Source)) + 
      scale_color_viridis(discrete=TRUE) + 
      theme_bw() +
      theme(axis.title.y=element_blank()) +
      labs(x = xlab, title = title) +
      geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
      #geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
      geom_point(aes(shape = Source)) +
      scale_shape_manual(values=1:nlevels(df$Source)) +
      theme(legend.position = 'None') +
      guides(color=guide_legend(override.aes=list(linetype=c(rep(0,1),rep(1,x$n_sources))))) +
      ggargs
  } else {
    g = ggplot(data=df, aes(x = x, y = y, colour = Source)) + 
      scale_color_grey() + 
      theme_bw() +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      labs(x = xlab, title = title) +
      geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
      #geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
      geom_point(aes(shape = Source)) +
      scale_shape_manual(values=1:nlevels(df$Source)) +
      theme(legend.title=element_blank(),legend.key = element_blank()) +
      guides(color=guide_legend(override.aes=list(linetype=c(rep(0,1),rep(1,x$n_sources))))) +
      ggargs
  }
}

print(g)

}
