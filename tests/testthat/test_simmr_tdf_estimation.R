set.seed(123)
library(vdiffr)
co <- function(expr) capture.output(expr, file = "NUL")

mix <- matrix(c(
  -10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54,
  -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89,
  11.73, 10.89, 11.05, 12.3
), ncol = 2, nrow = 10)
colnames(mix) <- c("d13C", "d15N")
s_names <- c("Source A", "Source B", "Source C", "Source D")
s_means <- matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol = 2, nrow = 4)
s_sds <- matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol = 2, nrow = 4)
conc <- matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol = 2, nrow = 4)

# Load into simmr
simmr_tdf <- simmr_load(
  mixtures = mix,
  source_names = s_names,
  source_means = s_means,
  source_sds = s_sds,
  concentration_means = conc
)

# MCMC run
co(simmr_tdf_out <- simmr_mcmc_tdf(simmr_tdf,
  p = matrix(rep(
    1 / simmr_tdf$n_sources,
    simmr_tdf$n_sources
  ),
  ncol = simmr_tdf$n_sources,
  nrow = simmr_tdf$n_obs, byrow = TRUE
  ),
  mcmc_control = list(iter = 100, burn = 10, thin = 1, n.chain = 4)
))

# Now put these corrections back into the model and check the
# iso-space plots and dietary output
simmr_tdf_2 <- simmr_load(
  mixtures = mix,
  source_names = s_names,
  source_means = s_means,
  source_sds = s_sds,
  correction_means = simmr_tdf_out$c_mean_est,
  correction_sds = simmr_tdf_out$c_sd_est,
  concentration_means = conc
)

test_that("main tdf function works", {
  co(tdf1 <- summary(simmr_tdf_out, type = "diagnostics"))
  expect_true(is.list(tdf1))
  expect_identical(names(tdf1), c("gelman", "quantiles", "statistics", "correlations"))
})

test_that("Other summary tdf functions produce output", {
  expect_output(summary(simmr_tdf_out, type = "quantiles"))
  expect_output(summary(simmr_tdf_out, type = "statistics"))
  expect_output(summary(simmr_tdf_out, type = "correlations"))
})

test_that("tdf output can be re-used", {
  co(simmr_tdf_2_out <- simmr_mcmc(simmr_tdf_2,
    mcmc_control = list(iter = 100, burn = 10, thin = 1, n.chain = 4)
  ))
  # Plot with corrections now
  p <- plot(simmr_tdf_2)
  expect_doppelganger('tdf_corrected_1', p)
  co(s1 <- summary(simmr_tdf_2_out, type = "diagnostics"))
  expect_true(is.list(s1))
  expect_identical(names(s1), c("gelman", "quantiles", "statistics", "correlations"))
  p2 <- plot(simmr_tdf_2_out, type = "boxplot")
  expect_doppelganger('tdf_corrected_2', p)
  
})
