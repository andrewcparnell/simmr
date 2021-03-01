set.seed(123)
co <- function(expr) capture.output(expr, file = "NUL")

data("geese_data_day1")
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

test_that("basic simmr_input plot", {
  p <- plot(simmr_1)
  expect_doppelganger("simmr_input", p)
  p <- plot(simmr_1, colour = FALSE)
  expect_doppelganger("simmr_input_no_col", p)
  p <- plot(simmr_1, tracers = c(2, 1))
  expect_doppelganger("simmr_input_rev_tracers", p)
})

test_that("1D simmr plot", {
  simmr_2 <- with(
    geese_data_day1,
    simmr_load(
      mixtures = mixtures[1, , drop = FALSE],
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = concentration_means
    )
  )

  p <- plot(simmr_2)
  expect_doppelganger("simmr_input_1obs", p)
  p <- plot(simmr_2, colour = FALSE)
  expect_doppelganger("simmr_input_1obs_nocol", p)
})

test_that("Multi-groups plot", {
  data(geese_data)
  simmr_4 <- with(
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

  p <- plot(simmr_4)
  expect_doppelganger("simmr_input_groups", p)
  p <- plot(simmr_4, group = 5)
  expect_doppelganger("simmr_input_groups_specified", p)
})

test_that("Single iso plot", {
  data(geese_data)
  simmr_5 <- with(
    geese_data,
    simmr_load(
      mixtures = mixtures[, 1, drop = FALSE],
      source_names = source_names,
      source_means = source_means[, 1, drop = FALSE],
      source_sds = source_sds[, 1, drop = FALSE],
      correction_means = correction_means[, 1, drop = FALSE],
      correction_sds = correction_sds[, 1, drop = FALSE],
      concentration_means = concentration_means[, 1, drop = FALSE],
      group = groups
    )
  )

  p <- plot(simmr_5)
  expect_doppelganger("simmr_input_iso1", p)
  p <- plot(simmr_5, group = 5)
  expect_doppelganger("simmr_input_iso1_group5", p)
  p <- plot(simmr_5, group = 5, colour = FALSE)
  expect_doppelganger("simmr_input_iso1_group5_nocol", p)
})
