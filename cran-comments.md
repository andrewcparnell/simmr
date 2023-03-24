## Test environments
* local OS X install, R 4.2.3
* GitHub actions
* R-hub
* win-builder (devel and release)

## R CMD check results
There were no notes, errors or warnings.

On win-builder I got a warning about a possibly invalid DOI in the DESCRIPTION. But when I checked this the DOI worked perfectly. 

One of the R-hub runs failed (Windows Server 2022, R-devel, 64 bit) because it couldn't find JAGS but I think this is a missing dependency on R-hub rather than a problem with the package. 

## Previous CRAN package checks
There are no previous errors or warnings for the package on CRAN.

## Reverse dependencies
There are currently no reverse dependencies