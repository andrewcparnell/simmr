set.seed(123)
library(vdiffr)
co <- function(expr) capture.output(expr, file = NULL)

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
co(simmr_1_out <- simmr_mcmc(simmr_1,
  mcmc_control = list(iter = 100, burn = 10, thin = 1, n.chain = 4)
))

co(simmr_1ffvb_out <- simmr_ffvb(simmr_1,
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

data(geese_data)
simmr_2 <- with(
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
co(simmr_2_out <- simmr_mcmc(simmr_2,
  mcmc_control = list(iter = 100, burn = 10, thin = 1, n.chain = 4)
))
co(simmr_2ffvb_out <- simmr_ffvb(simmr_2,
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

test_that("prior viz for 1 group", {
  p1 <- prior_viz(simmr_1_out)
  expect_list(p1)
  expect_class(p1$plot, "ggplot")
  expect_class(p1$p_prior_sim, "matrix")
  # Change some options
  p1a <- prior_viz(simmr_1_out, plot = FALSE, include_posterior = FALSE, n_sims = 10)
  expect_matrix(p1a)
})

test_that("prior viz for multiple groups", {
  p2 <- prior_viz(simmr_2_out)
  expect_list(p2)
  expect_class(p2$plot, "ggplot")
  expect_class(p2$p_prior_sim, "matrix")
  p3 <- prior_viz(simmr_2_out, group = 2)
  expect_list(p3)
  expect_false(p2$p_prior_sim[1, 1] == p3$p_prior_sim[1, 1])
  # Change some options
  p4 <- prior_viz(simmr_2_out, group = 2, plot = TRUE, include_posterior = FALSE, n_sims = 10)
  expect_list(p4)
  expect_error(prior_viz(simmr_2_out, group = 1.5, n_sims = 10))
  expect_error(prior_viz(simmr_2_out, group = 12, n_sims = 10))
})

test_that("prior viz for ffvb 1 group", {
  p1 <- prior_viz(simmr_1ffvb_out)
  expect_list(p1)
  expect_class(p1$plot, "ggplot")
  expect_class(p1$p_prior_sim, "matrix")
  # Change some options
  p1a <- prior_viz(simmr_1ffvb_out, plot = FALSE, include_posterior = FALSE, n_sims = 10)
  expect_matrix(p1a)
})

test_that("prior viz for ffvb for multiple groups", {
  p2 <- prior_viz(simmr_2ffvb_out)
  expect_list(p2)
  expect_class(p2$plot, "ggplot")
  expect_class(p2$p_prior_sim, "matrix")
  p3 <- prior_viz(simmr_2ffvb_out, group = 2)
  expect_list(p3)
  expect_false(p2$p_prior_sim[1, 1] == p3$p_prior_sim[1, 1])
  # Change some options
  p4 <- prior_viz(simmr_2ffvb_out, group = 2, plot = TRUE, include_posterior = FALSE, n_sims = 10)
  expect_list(p4)
  expect_error(prior_viz(simmr_2ffvb_out, group = 1.5, n_sims = 10))
  expect_error(prior_viz(simmr_2ffvb_out, group = 12, n_sims = 10))
})


test_that("posterior predictive for 1 groups", {
  co(pp1 <- posterior_predictive(simmr_1_out))
  expect_true(is.list(pp1))
  expect_true(is.data.frame(pp1$table))
  expect_true(is.numeric(pp1$p))
  # Change some options
  co(pp2 <- posterior_predictive(simmr_1_out, prob = 0.7, plot_ppc = FALSE))
  expect_true(is.list(pp2))
  expect_true(is.data.frame(pp2$table))
  expect_true(is.numeric(pp2$p))
})

test_that("posterior predictive for multiple groups", {
  co(pp3 <- posterior_predictive(simmr_2_out, group = 1))
  expect_true(is.list(pp3))
  expect_true(is.data.frame(pp3$table))
  expect_true(is.numeric(pp3$p))
  co(pp4 <- posterior_predictive(simmr_2_out, group = 2))
  expect_false(pp3$table[1, 1] == pp4$table[1, 1])

  # Change some options
  co(pp5 <- posterior_predictive(simmr_2_out, group = 2, prob = 0.7, plot_ppc = FALSE))
  expect_true(is.list(pp5))
  expect_true(is.data.frame(pp5$table))
  expect_true(is.numeric(pp5$p))
})

test_that("simmr elicit function", {
  co(np1 <- simmr_elicit(
    n_sources = 4,
    proportion_means = c(0.5, 0.2, 0.2, 0.1),
    proportion_sds = c(0.08, 0.02, 0.01, 0.02),
    n_sims = 10
  ))
  expect_true(is.list(np1))
  expect_true(length(np1$mean) == 4)
  expect_true(length(np1$sd) == 4)
  expect_warning(co(simmr_elicit(
    n_sources = 4,
    proportion_means = c(0.5, 0.2, 0.2, 0.1),
    proportion_sds = c(1, 1, 1, 0.02),
    n_sims = 10
  )))
  expect_error(co(simmr_elicit(
    n_sources = 4,
    proportion_means = c(-0.5, 0.2, 0.2, 0.1),
    proportion_sds = c(1, 1, 1, 0.02),
    n_sims = 10
  )))
  expect_error(co(simmr_elicit(
    n_sources = 4,
    proportion_means = c(0.5, 0.2, 0.2, 0.1),
    proportion_sds = c(1, 1.5, 1, 0.02),
    n_sims = 10
  )))
})
