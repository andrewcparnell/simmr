[![cran version](https://www.r-pkg.org/badges/version/simmr)](https://cran.rstudio.com/web/packages/simmr) 
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/simmr?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/simmr?color=82b4e8)](https://github.com/metacran/cranlogs.app)
[![Codecov test coverage](https://codecov.io/gh/andrewcparnell/simmr/branch/master/graph/badge.svg)](https://codecov.io/gh/andrewcparnell/simmr?branch=master)
[![R build status](https://github.com/andrewcparnell/simmr/workflows/R-CMD-check/badge.svg)](https://github.com/andrewcparnell/simmr/actions)

<a href="http://andrewcparnell.github.io/simmr/"><img src="https://raw.githubusercontent.com/andrewcparnell/simmr/master/badge/simmr_badge.png" height="200" align="right" /></a>


**Welcome to simmr!**

[simmr](http://andrewcparnell.github.io/simmr/) is a Bayesian stable isotope mixing model implemented in R which also uses [JAGS](https://mcmc-jags.sourceforge.net). It is intended as a replacement to the [SIAR](https://github.com/AndrewLJackson/siar) package. 

If you want the official stable version of the package from CRAN then go to R and type:

```
install.packages('simmr')
```

You can then load the package and view either the quick start or the full user manuals with:

```
library(simmr)
vignette("simmr")
vignette("quick_start")
```

Alternatively you can install the latest development version of the package in R with the commands:

```
if(!require(remotes)) install.packages('remotes')
remotes::install_github('andrewcparnell/simmr')
```

You can then load the package with

```
library(simmr)
```

... and look at the user manuals via:

```
vignette('simmr')
vignette('quick_start')
```
