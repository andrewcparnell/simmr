\name{plot.simmr_input}
\alias{plot.simmr_input}
\title{
Plot the \code{simmr_input} data created from \code{simmr_load}
}
\description{
This function creates iso-space (AKA tracer-space or delta-space) plots. They are vital in determining whether the data are suitable for running in a SIMM. 
}
\usage{
\method{plot}{simmr_input}(x, 
tracers = c(1, 2), 
title = "Tracers plot", 
xlab = "d13C", 
ylab = "d15N", 
sigmas = 1, 
colour = TRUE,...)
}

\arguments{
  \item{x}{
An object created via the function \code{\link{simmr_load}}
}
  \item{tracers}{
The choice of tracers to plot. If there are more than two tracers, it is recommended to plot every pair of tracers to determine whether the mixtures lie in the mixing polygon defined by the sources
}
  \item{title}{
A title for the graph
}
  \item{xlab}{
The x-axis label. By default this is assumed to be delta-13C but can be made richer if required. See examples below.
}
  \item{ylab}{
The y-axis label. By default this is assumed to be delta-15N in per mil but can be changed as with the x-axis label
}
  \item{sigmas}{
The number of standard deviations to plot on the source values. Defaults to 1.
}
  \item{colour}{
If TRUE (default) creates a plot. If not, puts the plot in black and white
}
  \item{...}{
Not used
}

}
\details{
It is desirable to have the vast majority of the mixture observations to be inside the convex hull defined by the food sources. When there are more than two tracers (as in one of the examples below) it is recommended to plot all the different pairs of the food sources. See the vignette for further details of richer plots
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
See \code{\link{plot.simmr_output}} for plotting the output of a simmr run. See \code{\link{simmr_mcmc}} for running a simmr object once the iso-space is deemed acceptable.
}
\examples{
# A simple example with 10 observations, 4 food sources and 2 tracers
mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54, 
-10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89, 
11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10)
colnames(mix) = c('d13C','d15N')
s_names=c('Source A','Source B','Source C','Source D')
s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, nrow=4)
s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4)
c_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol=2, nrow=4)
c_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol=2, nrow=4)
conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4)

simmr_1 = simmr_load(mixtures=mix,
                     source_names=s_names,
                     source_means=s_means,
                     source_sds=s_sds,
                     correction_means=c_means,
                     correction_sds=c_sds,
                     concentration_means = conc)

# Plot
plot(simmr_1)

### A more complicated example with 30 obs, 3 tracers and 4 sources
mix = matrix(c(-11.67, -12.55, -13.18, -12.6, -11.77, -11.21, -11.45, 
               -12.73, -12.49, -10.6, -12.26, -12.48, -13.07, -12.67, -12.26, 
               -13.12, -10.83, -13.2, -12.24, -12.85, -11.65, -11.84, -13.26, 
               -12.56, -12.97, -12.18, -12.76, -11.53, -12.87, -12.49, 7.79, 
               7.85, 8.25, 9.06, 9.13, 8.56, 8.03, 7.74, 8.16, 8.43, 7.9, 8.32, 
               7.85, 8.14, 8.74, 9.17, 7.33, 8.06, 8.06, 8.03, 8.16, 7.24, 7.24, 
               8, 8.57, 7.98, 7.2, 8.13, 7.78, 8.21, 11.31, 10.92, 11.3, 11, 
               12.21, 11.52, 11.05, 11.05, 11.56, 11.78, 12.3, 10.87, 10.35, 
               11.66, 11.46, 11.55, 11.41, 12.01, 11.97, 11.5, 11.18, 11.49, 
               11.8, 11.63, 10.99, 12, 10.63, 11.27, 11.81, 12.25), ncol=3, nrow=30)
colnames(mix) = c('d13C','d15N','d34S')
s_names = c('Source A', 'Source B', 'Source C', 'Source D') 
s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96, 
                   10.35, 7.51, 10.31, 9), ncol=3, nrow=4)
s_sds = matrix(c(0.46, 0.39, 0.42, 0.48, 0.44, 0.37, 0.49, 0.47, 0.49, 
                 0.42, 0.41, 0.42), ncol=3, nrow=4)
c_means = matrix(c(1.3, 1.58, 0.81, 1.7, 1.73, 1.83, 1.69, 3.2, 0.67, 
                   2.99, 3.38, 1.31), ncol=3, nrow=4)
c_sds = matrix(c(0.32, 0.64, 0.58, 0.46, 0.61, 0.55, 0.47, 0.45, 0.34, 
                 0.45, 0.37, 0.49), ncol=3, nrow=4)
conc = matrix(c(0.05, 0.1, 0.06, 0.07, 0.07, 0.03, 0.07, 0.05, 0.1, 
                0.05, 0.12, 0.11), ncol=3, nrow=4)

# Load this in:
simmr_3 = simmr_load(mixtures=mix,
                     source_names=s_names,
                     source_means=s_means,
                     source_sds=s_sds,
                     correction_means=c_means,
                     correction_sds=c_sds,
                     concentration_means = conc)


# Plot 3 times - first default d13C vs d15N 
plot(simmr_3)
# Now plot d15N vs d34S
plot(simmr_3,tracers=c(2,3))
# and finally d13C vs d34S
plot(simmr_3,tracers=c(1,3))
# See vignette('simmr') for fancier x-axis labels
}
