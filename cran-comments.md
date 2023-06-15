## Test environments
* local OS X install, R 4.2.3
* GitHub actions
* win-builder (devel and release)

## R CMD check results
There were no notes, errors or warnings.

On win-builder I got a warning about a possibly invalid DOI in the DESCRIPTION. But when I checked this the DOI worked perfectly. 

The R-hub runs failed because it couldn't find JAGS but I think this is a missing dependency on R-hub rather than a problem with the package (or a missing option that I need to add in somehow).

When I previously submitted the package there was a note about requiring C++11 which I have now fixed by removing that specification.

## Previous CRAN package checks
() added to all function names in the description texts
Added value to .Rd files to explain function results in the documentation.
Changed dontrun to donttest. (left as dontrun in compare_sources, 
plot.simmr_output prior_viz, simmr_elicit, summary.simmr_output, simmr_mcmc_tdf otherwise an error occurs with model_connection being 
left open - error with JAGS, and left as dontrun in simmr_mcmc and simmr_ffvb due to long times otherwise).
Changed cat/print to message.
Removed setting options and used supressWarnings instead

Package size was too large.
Version number was not updated.

## Reverse dependencies
There are currently no reverse dependencies