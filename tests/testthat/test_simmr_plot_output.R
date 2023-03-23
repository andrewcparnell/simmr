set.seed(123)
library(vdiffr)
co <- function(expr) capture.output(expr, file = "NUL")

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
set.seed(123)
co(simmr_out <- simmr_mcmc(simmr_in,
  mcmc_control = list(
    iter = 100,
    burn = 10,
    thin = 1,
    n.chain = 2
  )
))
set.seed(123)
co(simmr_out_ffvb <- simmr_ffvb(simmr_in,
  ffvb_control = list(
    n_output = 3600,
    S = 10,
    P = 1,
    beta_1 = 0.9,
    beta_2 = 0.9,
    tau = 1000,
    eps_0 = 0.1,
    t_W = 1
  )
))

# Taken from the simmr_mcmc example
test_that("plot.simmr_output", {
  p <- plot(simmr_out, type = "density")
  vdiffr::expect_doppelganger("plot_output_dens", p)
  p <- plot(simmr_out, type = "boxplot")
  vdiffr::expect_doppelganger("plot_output_box", p)
  p <- plot(simmr_out, type = "isospace")
  vdiffr::expect_doppelganger("out_iso", p)
  p <- plot(simmr_out, type = "histogram")
  vdiffr::expect_doppelganger("plot_output_hist", p)
  p <- plot(simmr_out, type = "matrix")
  vdiffr::expect_doppelganger("plot_output_matrix", p)
})
#
test_that("plot.simmr_output", {
  p <- plot(simmr_out_ffvb, type = "density")
  vdiffr::expect_doppelganger("plot_output_dens_ffvb", p)
  p <- plot(simmr_out_ffvb, type = "boxplot")
  vdiffr::expect_doppelganger("plot_output_box_ffvb", p)
  p <- plot(simmr_out, type = "isospace")
  vdiffr::expect_doppelganger("out_iso_ffvb", p)
  p <- plot(simmr_out_ffvb, type = "histogram")
  vdiffr::expect_doppelganger("plot_output_hist_ffvb", p)
  p <- plot(simmr_out_ffvb, type = "matrix")
  vdiffr::expect_doppelganger("plot_output_matrix_ffvb", p)
})
