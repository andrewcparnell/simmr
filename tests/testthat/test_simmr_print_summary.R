set.seed(123)
library(vdiffr)
co <- function(expr) capture.output(expr, file = "NUL")

data(geese_data_day1)

# Load into simmr
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
# MCMC run
co(simmr_1_out <- simmr_mcmc(simmr_1,
  mcmc_control = list(
    iter = 100,
    burn = 10,
    thin = 1,
    n.chain = 2
  )
))

# FFVB run
co(simmr_1_out_ffvb <-simmr_ffvb(simmr_1))

# Taken from the simmr_mcmc example
test_that("print.simmr_input", {
  expect_output(print(simmr_1))
})



test_that("print.simmr_output", {
  expect_output(print(simmr_1_out))
})

test_that("print.simmr_output", {
  expect_output(print(simmr_1_out_ffvb))
})

test_that("summary.simmr_output", {
  expect_output(summary(simmr_1_out, type = "diagnostics"))
  expect_output(summary(simmr_1_out))
})

test_that("summary.simmr_output", {
  expect_output(summary(simmr_1_out_ffvb))
  expect_output(summary(simmr_1_out_ffvb))
})


# Group version -----------------------------------------------------------

data(geese_data)
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
co(simmr_groups_out <- simmr_mcmc(simmr_groups,
  mcmc_control = list(
    iter = 100,
    burn = 10,
    thin = 1,
    n.chain = 2
  )
))
# Taken from the simmr_mcmc example
test_that("print.simmr_input", {
  expect_output(print(simmr_groups))
})

test_that("print.simmr_output", {
  expect_output(print(simmr_groups_out))
})

test_that("summary.simmr_output", {
  expect_output(summary(simmr_groups_out, type = "diagnostics"))
  expect_output(summary(simmr_groups_out))
})

test_that("summary.simmr_output groups ", {
  expect_output(summary(simmr_groups_out, type = "statistics", group = 1))
  expect_output(summary(simmr_groups_out, type = "statistics", group = c(1, 3)))
  co(summ_1 <- summary(simmr_groups_out, type = "statistics", group = 1))
  co(summ_2 <- summary(simmr_groups_out, type = "statistics", group = 2))
  expect_false(summ_1$statistics[[1]][1, 1] == summ_2$statistics[[1]][1, 1])
})
