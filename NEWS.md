# simmr 0.4.3.9000

  - Updated the simmr_elicit function to provide a more explicit warning for bad input objects
  - Updated compare_sources and compare_groups to allow exporting of the plot object for editing purposes

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

