\name{compare_sources}
\alias{compare_sources}
\title{
Compare dietary proportions between multiple sources
}
\description{
This function takes in an object of class \code{simmr_output} and creates probabilistic comparisons between the suppled sources. The group number can also be specified. 
}
\usage{
compare_sources(simmr_out, 
source_names=simmr_out$input$source_names,
group=1,
plot=TRUE)
}
\arguments{
  \item{simmr_out}{
An object of class \code{simmr_output} created from \code{\link{simmr_mcmc}}.
}
  \item{source_names}{
The names of at least two sources. These should match the names exactly given to \code{\link{simmr_load}}.
}
  \item{group}{
The integer values of the group numbers to be compared. If not specified assumes the first or only group
}
  \item{plot}{
A logical value specifying whether plots should be produced or not.
}
}
\details{
When two sources are specified, the function produces a direct calculation of the probability that the dietary proportion for one source is bigger than the other. When more than two sources are given, the function produces a set of most likely probabilistic orderings for each combination of sources. The function produces boxplots by default and also allows for the storage of the output for futher analysis if required.
}
\value{
If there are two sources, a vector containing the differences between the two dietary proportion proportions for these two sources. If there are multiple sources, a list containing the following fields:
\item{Ordering }{The different possible orderings of the dietary proportions across sources}
\item{out_all }{The dietary proportions for these sources specified as columns in a matrix}
}

\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
See \code{\link{simmr_mcmc}} for complete examples.
}
\examples{\dontrun{
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

# Compare two sources
compare_sources(simmr_1_out,sources=c('Zostera','U.lactuca'))

# Compare multiple sources
compare_sources(simmr_1_out)
}
}
