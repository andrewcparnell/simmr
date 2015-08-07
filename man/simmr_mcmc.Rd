\name{simmr_mcmc}
\alias{simmr_mcmc}
\alias{print.simmr_output}
\title{
Run a \code{simmr_input} object through the main simmr Markov chain Monte Carlo (MCMC) function
}
\description{
This is the main function of simmr. It takes a \code{simmr_input} object created via \code{\link{simmr_load}}, runs an MCMC to determine the dietary proportions, and then outputs a \code{simmr_output} object for further analysis and plotting via \code{\link{summary.simmr_output}} and \code{\link{plot.simmr_output}}.
}
\usage{
simmr_mcmc(simmr_in, 
prior.control=list(means=rep(0,simmr_in$n_sources), sd=rep(1,simmr_in$n_sources)), 
mcmc.control=list(iter=10000,burn=1000,thin=10,n.chain=4))
}
\arguments{
  \item{simmr_in}{
An object created via the function \code{\link{simmr_load}}
}
  \item{prior.control}{
A list of values including arguments named \code{means} and \code{sd} which represent the prior means and standard deviations of the dietary proportions in centralised log-ratio space. These can usually be left at their default values unless you wish to include to include prior information, in which case you should use the function \code{\link{simmr_elicit}}.
}
  \item{mcmc.control}{
A list of values including arguments named \code{iter} (number of iterations), \code{burn} (size of burn-in), \code{thin} (amount of thinning), and \code{n.chain} (number of MCMC chains). 
}
}
\details{
If, after running \code{\link{simmr_mcmc}} the convergence diagnostics in \code{\link{summary.simmr_output}} are not satisfactory, the values of \code{iter}, \code{burn} and \code{thin} in \code{mcmc.control} should be increased by a factor of 10.
}
\value{
An object of class \code{simmr_output} with two named top-level components:
\item{input }{The \code{simmr_input} object given to the \code{simmr_mcmc} function}
\item{output }{A set of MCMC chains of class \code{mcmc.list} from the coda package. These can be analysed using the \code{\link{summary.simmr_output}} and \code{\link{plot.simmr_output}} functions.}
}
\references{
Andrew C. Parnell, Donald L. Phillips, Stuart Bearhop, Brice X. Semmens, Eric J. Ward, Jonathan W. Moore, Andrew L. Jackson, Jonathan Grey, David J. Kelly, and Richard Inger. Bayesian stable isotope mixing models. Environmetrics, 24(6):387–399, 2013.

Andrew C Parnell, Richard Inger, Stuart Bearhop, and Andrew L Jackson. Source partitioning using stable isotopes: coping with too much variation. PLoS ONE, 5(3):5, 2010.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
\code{\link{simmr_load}} for creating objects suitable for this function, \code{\link{plot.simmr_input}} for creating isospace plots, \code{\link{summary.simmr_output}} for summarising output, and \code{\link{plot.simmr_output}} for plotting output.
}
\examples{\dontrun{
## See the package vignette for a detailed run through of these 4 examples

# Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep

# The data
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

# Load into simmr
simmr_1 = simmr_load(mixtures=mix,
                     source_names=s_names,
                     source_means=s_means,
                     source_sds=s_sds,
                     correction_means=c_means,
                     correction_sds=c_sds,
                     concentration_means = conc)

# Plot
plot(simmr_1)

# Print
simmr_1

# MCMC run
simmr_1_out = simmr_mcmc(simmr_1)

# Print it
print(simmr_1_out)

# Summary
summary(simmr_1_out)
summary(simmr_1_out,type='diagnostics')
summary(simmr_1_out,type='correlations')
summary(simmr_1_out,type='statistics')
ans = summary(simmr_1_out,type=c('quantiles','statistics'))

# Plot
plot(simmr_1_out)
plot(simmr_1_out,type='boxplot')
plot(simmr_1_out,type='histogram')
plot(simmr_1_out,type='density')
plot(simmr_1_out,type='matrix')

#####################################################################################

# A version with just one observation
simmr_2 = simmr_load(mixtures=mix[1,,drop=FALSE], # drop required to keep the mixtures as a matrix 
                     source_names=s_names,
                     source_means=s_means,
                     source_sds=s_sds,
                     correction_means=c_means,
                     correction_sds=c_sds,
                     concentration_means = conc)

# Plot
plot(simmr_2)

# MCMC run - automatically detects the single observation
simmr_2_out = simmr_mcmc(simmr_2)

# Print it
print(simmr_2_out)

# Summary
summary(simmr_2_out)
summary(simmr_2_out,type='diagnostics')
ans = summary(simmr_2_out,type=c('quantiles'))

# Plot
plot(simmr_2_out)
plot(simmr_2_out,type='boxplot')
plot(simmr_2_out,type='histogram')
plot(simmr_2_out,type='density')
plot(simmr_2_out,type='matrix')

#####################################################################################

# Data set 2: 3 isotopes (d13C, d15N and d34S), 30 observations, 4 sources

# The data
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

# Load into simmr
simmr_3 = simmr_load(mixtures=mix,
                     source_names=s_names,
                     source_means=s_means,
                     source_sds=s_sds,
                     correction_means=c_means,
                     correction_sds=c_sds,
                     concentration_means = conc)

# Get summary
print(simmr_3)

# Plot 3 times
plot(simmr_3)
plot(simmr_3,tracers=c(2,3))
plot(simmr_3,tracers=c(1,3))
# See vignette('simmr') for fancier axis labels

# MCMC run
simmr_3_out = simmr_mcmc(simmr_3)

# Print it
print(simmr_3_out)

# Summary
summary(simmr_3_out)
summary(simmr_3_out,type='diagnostics')
summary(simmr_3_out,type='quantiles')
summary(simmr_3_out,type='correlations')

# Plot
plot(simmr_3_out)
plot(simmr_3_out,type='boxplot')
plot(simmr_3_out,type='histogram')
plot(simmr_3_out,type='density')
plot(simmr_3_out,type='matrix')

#####################################################################################

# Data set 4 - identified by Fry (2014) as a failing of SIMMs
# See the vignette for more interpreation of these data and the output

# The data
mix = matrix(c(-14.65, -16.39, -14.5, -15.33, -15.76, -15.15, -15.73, 
               -15.52, -15.44, -16.19, 8.45, 8.08, 7.39, 8.68, 8.23, 7.84, 8.48, 
               8.47, 8.44, 8.37), ncol=2, nrow=10)
s_names = c('Source A', 'Source B', 'Source C', 'Source D') 
s_means = matrix(c(-25, -25, -5, -5, 4, 12, 12, 4), ncol=2, nrow=4)
s_sds = matrix(c(1, 1, 1, 1, 1, 1, 1, 1), ncol=2, nrow=4)

# Load into simmr - note no corrections or concentrations
simmr_4 = simmr_load(mixtures=mix,
                     source_names=s_names,
                     source_means=s_means,
                     source_sds=s_sds)

# Get summary
print(simmr_4)

# Plot 
plot(simmr_4)

# MCMC run - needs slightly longer
simmr_4_out = simmr_mcmc(simmr_4,
mcmc.control=list(iter=100000,burn=10000,thin=100,n.chain=4))

# Print it
print(simmr_4_out)

# Summary
summary(simmr_4_out)
summary(simmr_4_out,type='diagnostics')
ans = summary(simmr_4_out,type=c('quantiles','statistics'))

# Plot
plot(simmr_4_out)
plot(simmr_4_out,type='boxplot')
plot(simmr_4_out,type='histogram')
plot(simmr_4_out,type='density') # Look at the massive correlations here
plot(simmr_4_out,type='matrix')
}
}
