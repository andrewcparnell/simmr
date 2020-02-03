Updates for 0.4.2

- Added in a function simmr_groups_re which allows for a hierarchical model based on borrowing strengths between groups
- Updated all plotting summary and print functions to work with the above

Updates for 0.4.1

- Fixed some major bugs to plot.simmr_output and compare_groups which caused the wrong groups to be selected
- Fixed a minor dependency bug as no longer using coda

Updates for 0.4.0

- Included `prior_viz` to visualise and contrast the prior and posterior distributions
- Included `posterior_posterior_predictive` to visualise model fit using `bayesplot`
- Added `simmr_mcmc_tdf` to estimate trophic discrimination factors for known diet studies
- Updated `simmr_mcmc` to use R2jags 

