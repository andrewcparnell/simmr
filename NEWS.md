# simmr 0.5.1

  - Fixed bug introduced in 0.5.0 affecting the prior distribution on the residual standard deviation
  - Added extra argument to simmr_mcmc and simmr_mcmc_tdf to enable more control on the prior distribution on the residual standard deviation

# simmr 0.5.0

  - Implemented new Fixed Form Variational Bayes method for fitting SIMMs (publication forthcoming)
  - Added new GGally matrix plots

# simmr 0.4.6

  - Implemented vdiffr for better checking of output plots
  - Added capture.output to remove verbose testing
  - Added feature to allow for >2 sources to combined in combine_sources
  - Added ability to use data in matrix, data frame, or tibble format without error
  - Changed the way prior_viz object is plotted and returned to allow for greater customisation
  - Improved test that to test for bad source mean and sd shapes
  - Added a new vignette on advanced plotting (and moved other parts out of main vignette)

# simmr 0.4.5

  - Updated new checkmate error checking for multiple functions
  - Added new tests for 90%+ code coverage
  - Fixed a bug that stopped some plots being outputted correctly
  - Fixed some broken examples

# simmr 0.4.4

  - Updated the simmr_elicit function to provide a more explicit warning for bad input objects
  - Updated compare_sources and compare_groups to allow exporting of the plot object for editing purposes
  - Used styler to correct code style

# simmr 0.4.3

  - Fixed a bug in combine_sources that stopped it working for multiple groups
  - Added in a load more tests to increase code coverage
  - Implemented checkmate for error checking in simmr_load

# simmr 0.4.2

  - Fixed a bug with the summary function which was always reporting the same group when an individual group was specified (didn't apply to multiple group summaries). Added a test for that bug.
  - Updated posterior_predictive to produce some more helpful output

# simmr 0.4.1

  - Fixed some major bugs to plot.simmr_output and compare_groups which caused the wrong groups to be selected
  - Fixed a minor dependency bug as no longer using coda

# simmr 0.4.0

  - Included `prior_viz` to visualise and contrast the prior and posterior distributions
  - Included `posterior_posterior_predictive` to visualise model fit using `bayesplot`
  - Added `simmr_mcmc_tdf` to estimate trophic discrimination factors for known diet studies
  - Updated `simmr_mcmc` to use R2jags 

