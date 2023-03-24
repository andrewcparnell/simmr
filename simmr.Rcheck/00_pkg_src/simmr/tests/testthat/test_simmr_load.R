context("simmr_load")

library(simmr)
co <- function(expr) capture.output(expr, file = NULL)

test_that("simplest example", {
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
})

test_that("group example", {
  data("geese_data")
  simmr_groups <- with(
    geese_data,
    simmr_load(
      mixtures = mixtures,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = concentration_means,
      group = groups
    )
  )

  expect_s3_class(simmr_groups, "simmr_input")
  expect_true(is.matrix(simmr_groups$source_means))
  expect_true(is.matrix(simmr_groups$source_sds))
  expect_true(is.matrix(simmr_groups$correction_means))
  expect_true(is.matrix(simmr_groups$correction_sds))
  expect_true(is.matrix(simmr_groups$concentration_means))
  expect_true(is.factor(simmr_groups$group))
  expect_true(is.integer(simmr_groups$group_int))
})

# Would be good to test for a load of errors at this point

test_that("test without tefs and concentration", {
  data("geese_data_day1")

  # Load in with simmr_load
  simmr_2 <- with(
    geese_data_day1,
    simmr_load(
      mixtures = mixtures,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds
    )
  )
  expect_s3_class(simmr_2, "simmr_input")
  expect_true(is.matrix(simmr_2$source_means))
  expect_true(is.matrix(simmr_2$source_sds))
  expect_true(is.matrix(simmr_2$correction_means))
  expect_true(is.matrix(simmr_2$correction_sds))
  expect_true(is.matrix(simmr_2$concentration_means))
})

# Test if you're accidentally giving it text data

test_that("test without tefs and concentration", {
  data("geese_data_day1")
  geese_data_day1_tmp <- geese_data_day1
  geese_data_day1_tmp$mixtures <- matrix(as.character(geese_data_day1$mixtures), nrow = 9, ncol = 2)

  # Load in with simmr_load
  expect_error(with(
    geese_data_day1_tmp,
    simmr_load(
      mixtures = mixtures,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds
    )
  ))
})

test_that("test it works from excel", {
  library(readxl)
  path <- system.file("extdata", "geese_data.xls", package = "simmr")
  geese_data <- lapply(excel_sheets(path), read_excel, path = path)

  targets <- geese_data[[1]]
  sources <- geese_data[[2]]
  TEFs <- geese_data[[3]]
  concdep <- geese_data[[4]]

  geese_simmr <- simmr_load(
    mixtures = targets[, 1:2],
    source_names = sources$Sources,
    source_means = sources[, 2:3],
    source_sds = sources[, 4:5],
    correction_means = TEFs[, 2:3],
    correction_sds = TEFs[, 4:5],
    concentration_means = concdep[, 2:3],
    group = as.factor(paste("Day", targets$Time))
  )
  co(simmr_out <- simmr_mcmc(geese_simmr,
    mcmc_control = list(
      iter = 100,
      burn = 10,
      thin = 1,
      n.chain = 2
    )
  ))
  expect_s3_class(simmr_out, "simmr_output")
})
