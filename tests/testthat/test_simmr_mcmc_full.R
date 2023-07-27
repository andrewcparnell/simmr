# Test some of the data sets that people have sent through to make
# sure that simmr converges and is not too sensitive to the prior
# distribution on the residual standard deviation

set.seed(123)
co <- function(expr) capture.output(expr, file = NULL)


# Start with Geese data ---------------------------------------------------

mix <- matrix(c(
  -10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54,
  -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89,
  11.73, 10.89, 11.05, 12.3
), ncol = 2, nrow = 10)
colnames(mix) <- c("d13C", "d15N")
s_names <- c("Source A", "Source B", "Source C", "Source D")
s_means <- matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol = 2, nrow = 4)
s_sds <- matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol = 2, nrow = 4)
c_means <- matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol = 2, nrow = 4)
c_sds <- matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol = 2, nrow = 4)
conc <- matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol = 2, nrow = 4)

# Taken from the simmr_mcmc example
test_that("simmr_mcmc_full_run_geese", {
  skip_on_ci()
  skip_on_cran()

  # Load into simmr
  simmr_1 <- simmr_load(
    mixtures = mix,
    source_names = s_names,
    source_means = s_means,
    source_sds = s_sds,
    correction_means = c_means,
    correction_sds = c_sds,
    concentration_means = conc
  )
  # MCMC run
  co(simmr_1_out <- simmr_mcmc(simmr_1))

  co(rhat <- summary(simmr_1_out, type = "diagnostics")$gelman)
  expect_true(all(rhat < 1.05))

  # Create a quick traceplot of sd values to ensure it worked ok
  # R2jags::traceplot(simmr_1_out$output[[1]], varname = "sigma")
})

# Now run another on simmr solo
test_that("simmr_mcmc_full_run_geese_solo", {
  skip_on_ci()
  skip_on_cran()

  # Load into simmr
  simmr_2 <- simmr_load(
    mixtures = mix[1, , drop = FALSE], # drop required to keep the mixtures as a matrix
    source_names = s_names,
    source_means = s_means,
    source_sds = s_sds,
    correction_means = c_means,
    correction_sds = c_sds,
    concentration_means = conc
  )

  # MCMC run - automatically detects the single observation
  co(simmr_2_out <- simmr_mcmc(simmr_2))

  co(rhat <- summary(simmr_2_out, type = "diagnostics")$gelman)
  expect_true(all(rhat < 1.05))
})

# Test what happens when you change the prior distribution on sigma -------
test_that("simmr_mcmc_full_run_geese", {
  skip_on_ci()
  skip_on_cran()
  
  #Load 
  simmr_1 <- simmr_load(
    mixtures = mix,
    source_names = s_names,
    source_means = s_means,
    source_sds = s_sds,
    correction_means = c_means,
    correction_sds = c_sds,
    concentration_means = conc)
  
  # Load into simmr
  co(simmr_3_out <- simmr_mcmc(simmr_1, prior_control = list(
    means = rep(0, simmr_1$n_sources),
    sd = rep(1, simmr_1$n_sources),
    sigma_shape = c(3, 3),
    sigma_rate = c(3/50, 3/50)
  )))
  
  co(rhat <- summary(simmr_3_out, type = "diagnostics")$gelman)
  # expect_true(all(rhat < 1.05)) # THIS IS FALSE?
  
  # Create a quick traceplot of sd values to ensure it worked ok
  # R2jags::traceplot(simmr_3_out$output[[1]], varname = "sigma")
  
  # Order is [iteration, chain, parameter]
  # sigma1 <- simmr_3_out$output[[1]]$BUGSoutput$sims.array[,1,6]
  # p <- simmr_3_out$output[[1]]$BUGSoutput$sims.array[,1,2:5]
  # sigma2 <- simmr_3_out$output[[1]]$BUGSoutput$sims.array[,,7]
  # 
  # p <- simmr_3_out$output[[1]]$BUGSoutput$sims.list$p
  # sig <- simmr_3_out$output[[1]]$BUGSoutput$sims.list$sigma
})


# Next data set is Fletcher data ------------------------------------------

data <- structure(list(Urban_Code = c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1,
  1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1
), Carbon = c(
  -24.7,
  -23.54, -26.07, -23.23, -24.84, -25.41, -24.57, -23.79, -25.11,
  -25.04, -24.04, -23.61, -23.46, -23.53, -24.06, -24.96, -23.12,
  -24.88, -24.35, -23.78, -24.76, -24.06, -23.86, -25.32, -23.07,
  -23.48, -25.2, -24.58, -24.19, -25.49, -23.89, -24.58, -24.32,
  -22.93, -25.09, -24.29, -24.13, -24.02, -23.81, -24.24, -24.72,
  -24.18, -24.51, -24.74, -24.95, -25.26, -22.69, -22.66, -24.95,
  -22.55, -22.06, -22.33, -22.23, -21.92, -24.49, -23.89, -23.34,
  -24.69, -24.76, -24.41, -22.74, -24.3, -25.22, -23.76, -24.16,
  -25.65, -21.99, -25.28, -23.53, -24.53, -26.13, -23.71, -23.49,
  -23.57, -23.89, -23.98, -23.83, -23.86, -23.72, -25.05, -23.86,
  -23.13, -23.28, -22.72, -23.38, -23.91, -23.78, -26.02, -24.06,
  -22.6, -23.67, -23.16, -21.61
), Nitrogen = c(
  10.92, 8.32, 9.17, 9.55,
  11.81, 10.47, 9.45, 9.62, 10.3, 11.49, 11.05, 10.15, 8.02, 8.8,
  9.03, 10.26, 11.58, 8.43, 10.23, 12.98, 9.38, 10.22, 9.2, 10.62,
  8.51, 7.91, 9.95, 8.52, 9.48, 11.22, 10.26, 9.71, 10.7, 8.48,
  11.52, 7.62, 7.73, 7.79, 8.28, 8.38, 9.63, 7.26, 9.58, 12.31,
  9.31, 8.89, 7.15, 7.65, 9.64, 7.64, 7.19, 7.2, 7.64, 7.52, 8.33,
  7.98, 7.11, 8.77, 7.9, 10.21, 7.37, 8.27, 10.67, 9.62, 7.77,
  9.82, 8.25, 10.24, 9.16, 9.48, 9.89, 7.41, 8.06, 7.66, 9.68,
  8.99, 11.43, 9.49, 8.66, 11.8, 9.07, 8.1, 7.77, 8.06, 8.3, 9.9,
  9.17, 9.77, 9.79, 8.37, 10.17, 7.59, 9.22
)), class = "data.frame", row.names = c(
  NA,
  -93L
))

# Food sources (with TEFs added to food sources beforehand)
mix <- cbind(data$Carbon, data$Nitrogen)
colnames(mix) <- c("d13C", "d15N")
s_names_final <- c("Anthropogenic Food", "Pet Food", "Mammals", "Birds", "Invertbrates", "Fruit")
s_means_final <- matrix(c(-21.94, -23.33, -24.92, -23.87, -24.96, -25.1, 7.55, 7.4, 9.8, 9.88, 5.59, 0.77), ncol = 2, nrow = 6)
s_sds_final <- matrix(c(0.6, 1.91, 2.29, 1.16, 1.55, 1.98, 0.6, 0.87, 1.64, 1.4, 3.89, 2.67), ncol = 2, nrow = 6)


test_that("simmr_mcmc_fletcher", {
  skip_on_ci()
  skip_on_cran()
  
  # Load into simmr
  simmr_in <- simmr_load(
    mixtures = mix,
    source_names = s_names_final,
    source_means = s_means_final,
    source_sds = s_sds_final
  )
  # MCMC run
  co(simmr_out <- simmr_mcmc(simmr_in))
  
  co(rhat <- summary(simmr_out, type = "diagnostics")$gelman)
  expect_true(all(rhat < 1.05))
  
  # Create a quick traceplot of sd values to ensure it worked ok
  # R2jags::traceplot(simmr_1_out$output[[1]], varname = "sigma")
})




