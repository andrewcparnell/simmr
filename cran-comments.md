## Test environments
* local OS X install, R 4.2.3
* GitHub actions
* win-builder (devel and release)

## R CMD check results
There were no notes, errors or warnings.

On win-builder I got a warning about a possibly invalid DOI in the DESCRIPTION. But when I checked this the DOI worked perfectly. 

The R-hub runs failed because it couldn't find JAGS but I think this is a missing dependency on R-hub rather than a problem with the package (or a missing option that I need to add in somehow).

## Previous CRAN package checks
There are no previous errors or warnings for the package on CRAN.

## Reverse dependencies
There are currently no reverse dependencies