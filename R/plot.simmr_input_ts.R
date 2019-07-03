#' Plot a simmr time series
#' 
#' This function creates iso-space (AKA tracer-space or delta-space) plots for time series models.
#' They are vital in determining whether the data are suitable for running in a
#' SIMM.
#' 
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
#' @seealso See \code{\link{plot.simmr_input}} for plotting standard simmr objects. See \code{\link{simmr_mcmc}} for running a time series simmr object once the iso-space is deemed acceptable.
#' @examples
#' 
#' @export
plot.simmr_input_ts <-
function(x,
         tracers = c(1,2),
         title = 'Tracers plot',
         xlab = x$tracer_names[tracers[1]],
         ylab = x$tracer_names[tracers[2]], 
         sigmas = 1, 
         mix_name = 'Mixtures', 
         time_units = 'Time',
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
size=c(rep(0.5,x$n_sources),x$time_var)

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
      geom_linerange(aes(x=x,ymax=y_upper,ymin=y_lower)) +
      geom_point(aes(x = x, y = y, shape = Source, size = size)) +
      scale_shape_manual(values=1:nlevels(df$Source)) +
      theme(legend.key = element_blank()) +
      labs(size = time_units) +
      guides(color=guide_legend(override.aes=list(linetype=c(rep(0,1),rep(1,x$n_sources))))) +
      ggargs
  } else {
    g=ggplot(data=df, aes(x = x,y = y,colour=Source)) + 
      theme_bw() +
      labs(x=xlab,y=ylab,title=title) +
      geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0)) +
      #geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source)) +
      geom_linerange(aes(x=x,ymax=y_upper,ymin=y_lower)) +
      geom_point(aes(x = x, y = y, shape = Source, size = size)) +
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
      geom_point(aes(shape = Source, size = size)) +
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
      geom_point(aes(shape = Source, size = size)) +
      scale_shape_manual(values=1:nlevels(df$Source)) +
      theme(legend.title=element_blank(),legend.key = element_blank()) +
      guides(color=guide_legend(override.aes=list(linetype=c(rep(0,1),rep(1,x$n_sources))))) +
      ggargs
  }
}

print(g)

}
