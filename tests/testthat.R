library(testthat)
library(simmr)
library(vdiffr)

co <- function(expr) capture.output(expr, file = "NUL")

test_check("simmr")
