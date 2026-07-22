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
  # Run the simmr_mcmc function with this informative prior
  co(simmr_1a_out <- simmr_mcmc(simmr_1,
                            prior_control=list(means = np1$mean,
                                               sd = np1$sd,
                                               sigma_shape = c(3,3), 
                                               sigma_rate = c(3/50, 3/50))))
  expect_list(simmr_1a_out)
  
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
  expect_error(co(simmr_elicit(
    n_sources = 4,
    proportion_means = c(0.5, 0.2, 0.2, 0.1),
    proportion_sds = c(0, 0, 0, 0),
    n_sims = 10
  )))
  expect_error(co(simmr_elicit(
    n_sources = 4,
    proportion_means = c(0.5, 0.2, 0.2, 0.1),
    proportion_sds = c(0, 0.02, 0.01, 0.02),
    n_sims = 10
  )))
})

# New posterior predictive test from Andrew J
all_simmr <- structure(list(mixtures = structure(list(`13C` = c(-18.2, -18.1532, 
                                                   -18.7622351, -17.9211375, -17.9917462, -18.1, -17.63575, -18.4474511, 
                                                   -18.0250179, -17.559484, -18.1533248, -18.5, -18.1402696, -18.1838125, 
                                                   -18.171362, -18.5443846, -18.7045555, -18.1603328, -18.8726632, 
                                                   -18.1, -18.1, -18.0578685, -18.0737192, -17.8743051, -18.0418088, 
                                                   -18.3596558, -17.9562178, -17.6052454, -18.523196, -18.3035493, 
                                                   -17.6704336, -17.8121304, -18.0961323, -18.5886016, -18.5005709, 
                                                   -18.4171348, -17.9496989, -18.1287323, -18.5089909, -18.519948, 
                                                   -18.3677551, -18.671112, -18.20455, -18.0634212, -18.8028064, 
                                                   -18.214656, -18.6442824, -18.5826214, -18.0270539, -18.7111375, 
                                                   -18.5888207, -18.745212, -18.1287, -19.1633275, -18.6413292, 
                                                   -18.3618337, -17.6443574, -18.4636087, -18.7771439, -18.3381481, 
                                                   -18.6117972, -18.4296832, -18.013458, -18.388544, -18.5897875, 
                                                   -18.5783167, -18.446418, -18.720512, -18.7998706, -18.7127752, 
                                                   -18.5376058, -19.2010966, -20.4047746, -18.931003, -18.5464132, 
                                                   -20.1845896), `15N` = c(8, 9.1232453, 8.2430793, 7.33993, 7.4226738, 
                                                                           7.4, 7.7564708, 7.370587, 8.1145779, 7.7013132, 8.1442693, 9.1, 
                                                                           8.01662, 9.0501504, 9.0022986, 9.1989609, 9.5876497, 7.9371352, 
                                                                           8.830732, 9, 7.9, 8.4560451, 8.7950948, 8.4309876, 8.3128379, 
                                                                           8.6780034, 8.8609743, 8.2295687, 9.0098449, 8.3881767, 7.7537447, 
                                                                           8.6404683, 8.6092366, 8.3733673, 8.7904867, 8.7024188, 7.5909333, 
                                                                           8.3320083, 6.34709, 10.615848, 10.7262546, 10.7940067, 10.5971589, 
                                                                           9.712436, 10.8625576, 9.6135809, 9.512452, 9.1437262, 9.656492, 
                                                                           9.6198901, 9.6658212, 10.7502135, 10.0027432, 8.085588, 9.781188, 
                                                                           9.0674481, 9.4352922, 8.9313896, 8.9957656, 8.739696, 8.88606, 
                                                                           9.105396, 8.7066201, 8.7715571, 9.4879944, 9.387964, 9.516404, 
                                                                           9.3667465, 9.6617437, 8.3079913, 8.9730271, 8.6040013, 8.7608866, 
                                                                           7.2620893, 10.558654, 8.979934)), row.names = c(NA, -76L), class = c("tbl_df", 
                                                                                                                                                "tbl", "data.frame")), source_names = c("Caribou", "Horse", "Mammoth", 
                                                                                                                                                                                        "Rodents", "Steppe-bison"), source_means = structure(list(mean13C = c(-18.89, 
                                                                                                                                                                                                                                                              -20.4811250916372, -20.7463507402576, -22.404347826087, -19.6802857142857
                                                                                                                                                                                        ), mean15N = c(3.85666666666667, 4.39463890659244, 6.81995522661881, 
                                                                                                                                                                                                       3.09565217391304, 4.78291964285714)), row.names = c(NA, -5L), class = c("tbl_df", 
                                                                                                                                                                                                                                                                               "tbl", "data.frame")), source_sds = structure(list(SD13C = c(0.534602656184946, 
                                                                                                                                                                                                                                                                                                                                            0.24513910866367, 0.275386736570222, 1.39690510540902, 0.302357672860575
                                                                                                                                                                                                                                                                               ), SD15N = c(1.4083134120879, 1.30845408820857, 0.565454516322188, 
                                                                                                                                                                                                                                                                                            1.4072792127395, 0.612842496277393)), row.names = c(NA, -5L), class = c("tbl_df", 
                                                                                                                                                                                                                                                                                                                                                                    "tbl", "data.frame")), correction_means = structure(list(mean13C = c(1.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                         1.1, 1.1, 1.1, 1.1), mean15N = c(3.8, 3.8, 3.8, 3.8, 3.8)), row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   -5L), class = c("tbl_df", "tbl", "data.frame")), correction_sds = structure(list(
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     SD13C = c(0.6, 0.6, 0.6, 0.6, 0.6), SD15N = c(0.8, 0.8, 0.8, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   0.8, 0.8)), row.names = c(NA, -5L), class = c("tbl_df", "tbl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 "data.frame")), concentration_means = structure(c(1, 1, 1, 1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   1, 1, 1, 1, 1, 1), dim = c(5L, 2L)), group = structure(c(1L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), levels = c("Lion", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "Short-faced bear", "Wolf"), class = "factor"), group_int = c(1L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), n_obs = 76L, n_tracers = 2L, 
               n_sources = 5L, n_groups = 3L), class = "simmr_input")


co(all_simmr_out <- simmr_mcmc(all_simmr))

test_that("posterior predictive for tibbles", {
  co(pp4 <- posterior_predictive(all_simmr_out, group = 1))
  expect_true(is.list(pp4))
  expect_true(is.data.frame(pp4$table))
  expect_true(is.numeric(pp4$p))
  expect_false(ncol(pp4$table) == 5)
})
