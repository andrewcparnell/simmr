Welcome to simmr!

simmr is a Bayesian stable isotope mixing model implemented in R which also uses [JAGS](http://mcmc-jags.sourceforge.net). It is intended as a replacement to the [SIAR](https://github.com/AndrewLJackson/siar) package. You can install the development version of the package in R with the commands:

```
library(devtools)
install_github('andrewcparnell/simmr',build_vignettes =TRUE)
```

You can then load the package with

```
library(simmr)
```

A vignette (i.e. a simple user manual) is available via:
```
vignette('simmr')


If you want the official stable version of the package from CRAN then go to R and type:

```
install.packages('simmr')
```

You can then load the package and view the vignette with:

```
library(simmr)
vignettes("simmr")
```
