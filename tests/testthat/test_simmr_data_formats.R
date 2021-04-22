context("simmr run on different data formats")

library(simmr)
library(tibble)

# Simmr should really work with data that are specified as matrices, data frames or tibbles. None of these should cause any problems

test_that("matrices", {
  data("geese_data_day1")

  # Load in with simmr_load
  simmr_1 <- with(
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
  expect_s3_class(simmr_1, "simmr_input")
  expect_true(is.matrix(simmr_1$source_means))
  expect_true(is.matrix(simmr_1$source_sds))
  expect_true(is.matrix(simmr_1$correction_means))
  expect_true(is.matrix(simmr_1$correction_sds))
  expect_true(is.matrix(simmr_1$concentration_means))
  expect_s3_class(simmr_mcmc(simmr_1), "simmr_output")
})

test_that("data frames", {
  data("geese_data_day1")
  
  # Load in with simmr_load
  simmr_1 <- with(
    geese_data_day1,
    simmr_load(
      mixtures = as.data.frame(mixtures),
      source_names = source_names,
      source_means = as.data.frame(source_means),
      source_sds = as.data.frame(source_sds),
      correction_means = as.data.frame(correction_means),
      correction_sds = as.data.frame(correction_sds),
      concentration_means = as.data.frame(concentration_means)
    )
  )
  expect_s3_class(simmr_1, "simmr_input")
  expect_true(is.data.frame(simmr_1$source_means))
  expect_true(is.data.frame(simmr_1$source_sds))
  expect_true(is.data.frame(simmr_1$correction_means))
  expect_true(is.data.frame(simmr_1$correction_sds))
  expect_true(is.data.frame(simmr_1$concentration_means))
  expect_s3_class(simmr_mcmc(simmr_1), "simmr_output")
})

test_that("tibbles", {
  data("geese_data_day1")
  
  # Load in with simmr_load
  simmr_1 <- with(
    geese_data_day1,
    simmr_load(
      mixtures = as_tibble(mixtures),
      source_names = source_names,
      source_means = as_tibble(source_means),
      source_sds = as_tibble(source_sds),
      correction_means = as_tibble(correction_means),
      correction_sds = as_tibble(correction_sds),
      concentration_means = as_tibble(concentration_means)
    )
  )
  expect_s3_class(simmr_1, "simmr_input")
  expect_true(is_tibble(simmr_1$source_means))
  expect_true(is_tibble(simmr_1$source_sds))
  expect_true(is_tibble(simmr_1$correction_means))
  expect_true(is_tibble(simmr_1$correction_sds))
  expect_true(is_tibble(simmr_1$concentration_means))
  expect_s3_class(simmr_mcmc(simmr_1), "simmr_output")
})
