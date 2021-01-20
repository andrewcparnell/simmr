context("simmr plot output")

library(simmr)

data(geese_data_day1)

# Load into simmr
simmr_in <- with(
  geese_data_day1,
  simmr_load(
    mixtures = mixtures,
    source_names = source_names,
    source_means = source_means,
    source_sds = source_sds,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means
  )
)

# MCMC run
simmr_out <- simmr_mcmc(simmr_in,
  mcmc_control = list(
    iter = 100,
    burn = 10,
    thin = 1,
    n.chain = 2
  )
)

# Taken from the simmr_mcmc example
test_that("plot.simmr_output", {
  expect_s3_class(plot(simmr_out, type = "density"), "ggplot")
  expect_s3_class(plot(simmr_out, type = "boxplot"), "ggplot")
  # expect_s3_class(plot(simmr_out,type='isospace'), 'ggplot')
  expect_s3_class(plot(simmr_out, type = "histogram"), "ggplot")
  expect_null(plot(simmr_out, type = "matrix"), "ggplot")
})
