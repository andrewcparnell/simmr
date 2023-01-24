set.seed(123)
co <- function(expr) capture.output(expr, file = "NUL")

# Single group version first
data("geese_data_day1")

# Taken from the simmr_mcmc example
test_that("simmr combine sources single group", {
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
  co(simmr_1_out_ffvb <- simmr_ffvb(simmr_1))

  # Combine two sources
  simmr_out_combine <- combine_sources(simmr_1_out,
    to_combine = c("Zostera", "Grass"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_combine, "simmr_output")
  expect_s3_class(simmr_out_combine$input, "simmr_input")
  expect_true(length(simmr_out_combine$input$source_names) == 3)
  expect_true(nrow(simmr_out_combine$input$correction_means) == 3)
  expect_true(nrow(simmr_out_combine$input$source_sds) == 3)
  expect_true(simmr_out_combine$input$n_sources == 3)
  expect_output(summary(simmr_out_combine))

  simmr_out_combine_ffvb <- combine_sources(simmr_1_out_ffvb,
    to_combine = c("Zostera", "Grass"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_combine_ffvb, "simmr_output")
  expect_s3_class(simmr_out_combine_ffvb$input, "simmr_input")
  expect_true(length(simmr_out_combine_ffvb$input$source_names) == 3)
  expect_true(nrow(simmr_out_combine_ffvb$input$correction_means) == 3)
  expect_true(nrow(simmr_out_combine_ffvb$input$source_sds) == 3)
  expect_true(simmr_out_combine_ffvb$input$n_sources == 3)
  expect_output(summary(simmr_out_combine_ffvb))

  # Check it works for multiple sources
  simmr_out_combine_mult <- combine_sources(simmr_1_out,
    to_combine = c("Zostera", "Grass", "Enteromorpha"),
    new_source_name = "Z, G, and E"
  )

  # Check it works for multiple sources
  simmr_out_combine_mult_ffvb <- combine_sources(simmr_1_out,
    to_combine = c("Zostera", "Grass", "Enteromorpha"),
    new_source_name = "Z, G, and E"
  )

  expect_s3_class(simmr_out_combine_mult, "simmr_output")
  expect_s3_class(simmr_out_combine_mult$input, "simmr_input")
  expect_true(length(simmr_out_combine_mult$input$source_names) == 2)
  expect_true(nrow(simmr_out_combine_mult$input$correction_means) == 2)
  expect_true(nrow(simmr_out_combine_mult$input$source_sds) == 2)
  expect_true(simmr_out_combine_mult$input$n_sources == 2)
  expect_output(summary(simmr_out_combine_mult))

  expect_s3_class(simmr_out_combine_mult_ffvb, "simmr_output")
  expect_s3_class(simmr_out_combine_mult_ffvb$input, "simmr_input")
  expect_true(length(simmr_out_combine_mult_ffvb$input$source_names) == 2)
  expect_true(nrow(simmr_out_combine_mult_ffvb$input$correction_means) == 2)
  expect_true(nrow(simmr_out_combine_mult_ffvb$input$source_sds) == 2)
  expect_true(simmr_out_combine_mult_ffvb$input$n_sources == 2)
  expect_output(summary(simmr_out_combine_mult_ffvb))

  # Haven't got an example with 5 sources to check?

  # Check some errors
  expect_error(combine_sources(simmr_1_out,
    to_combine = c("Zostera2", "Grass"),
    new_source_name = "Z and G"
  ))

  expect_error(combine_sources(simmr_1_out_ffvb,
    to_combine = c("Zostera2", "Grass"),
    new_source_name = "Z and G"
  ))

  expect_error(combine_sources(simmr_1_out,
    to_combine = c("Zostera", "Gass"),
    new_source_name = "Z and G"
  ))


  expect_error(combine_sources(simmr_1_out_ffvb,
    to_combine = c("Zostera", "Gass"),
    new_source_name = "Z and G"
  ))
})

# Taken from the simmr_mcmc example
test_that("simmr combine sources multiple group", {
  data("geese_data")

  # Load into simmr
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
  # MCMC run
  co(simmr_2_out <- simmr_mcmc(simmr_2,
    mcmc_control = list(
      iter = 100,
      burn = 10,
      thin = 1,
      n.chain = 2
    )
  ))

  # Combine two sources
  simmr_out_2_combine <- combine_sources(simmr_2_out,
    to_combine = c("Zostera", "Grass"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_2_combine, "simmr_output")
  expect_s3_class(simmr_out_2_combine$input, "simmr_input")
  expect_true(length(simmr_out_2_combine$input$source_names) == 3)
  expect_true(nrow(simmr_out_2_combine$input$correction_means) == 3)
  expect_true(nrow(simmr_out_2_combine$input$source_sds) == 3)
  expect_true(simmr_out_2_combine$input$n_sources == 3)
  expect_true(simmr_out_2_combine$input$n_groups == 8)

  # Check that it still works when you reverse the sources (previously crashed)
  simmr_out_3_combine <- combine_sources(simmr_2_out,
    to_combine = c("Grass", "Zostera"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_3_combine, "simmr_output")

  # Make sure the summaries are different
  co(summ_1 <- summary(simmr_out_2_combine, type = "statistics", group = 1))
  co(summ_2 <- summary(simmr_out_2_combine, type = "statistics", group = 2))
  expect_false(summ_1$statistics[[1]][1, 1] == summ_2$statistics[[1]][1, 1])

  # Make sure compare groups is different
  co(cg_1 <- compare_groups(simmr_out_2_combine, source_name = "Z and G", groups = 1:2))
  co(cg_2 <- compare_groups(simmr_out_2_combine, source_name = "Z and G", groups = 2:3))
  expect_false(mean(cg_1[[1]]) == mean(cg_2[[1]]))

  # Make sure compare sources is different
  co(cs_1 <- compare_sources(simmr_out_2_combine,
    source_names = c("Z and G", "Enteromorpha"),
    group = 1
  ))
  co(cs_2 <- compare_sources(simmr_out_2_combine,
    source_names = c("Z and G", "Enteromorpha"),
    group = 2
  ))
  expect_false(mean(cs_1[[1]]) == mean(cs_2[[1]]))

  # Make sure the plots are different
  cp_1 <- plot(simmr_out_2_combine, type = "histogram", group = 1)
  cp_2 <- plot(simmr_out_2_combine, type = "histogram", group = 2)
  expect_false(cp_1$data$Proportion[1] == cp_2$data$Proportion[1])
})

test_that("simmr combine sources multiple group", {
  data("geese_data")

  # Load into simmr
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
  # MCMC run
  co(simmr_2_out <- simmr_ffvb(simmr_2))

  # Combine two sources
  simmr_out_2_combine <- combine_sources(simmr_2_out,
    to_combine = c("Zostera", "Grass"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_2_combine, "simmr_output")
  expect_s3_class(simmr_out_2_combine$input, "simmr_input")
  expect_true(length(simmr_out_2_combine$input$source_names) == 3)
  expect_true(nrow(simmr_out_2_combine$input$correction_means) == 3)
  expect_true(nrow(simmr_out_2_combine$input$source_sds) == 3)
  expect_true(simmr_out_2_combine$input$n_sources == 3)
  expect_true(simmr_out_2_combine$input$n_groups == 8)

  # Check that it still works when you reverse the sources (previously crashed)
  simmr_out_3_combine <- combine_sources(simmr_2_out,
    to_combine = c("Grass", "Zostera"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_3_combine, "simmr_output")

  # Make sure the summaries are different
  co(summ_1 <- summary(simmr_out_2_combine, type = "statistics", group = 1))
  co(summ_2 <- summary(simmr_out_2_combine, type = "statistics", group = 2))
  expect_false(summ_1$statistics[[1]][1, 1] == summ_2$statistics[[1]][1, 1])

  # Make sure compare groups is different
  co(cg_1 <- compare_groups(simmr_out_2_combine, source_name = "Z and G", groups = 1:2))
  co(cg_2 <- compare_groups(simmr_out_2_combine, source_name = "Z and G", groups = 2:3))
  expect_false(mean(cg_1[[1]]) == mean(cg_2[[1]]))

  # Make sure compare sources is different
  co(cs_1 <- compare_sources(simmr_out_2_combine,
    source_names = c("Z and G", "Enteromorpha"),
    group = 1
  ))
  co(cs_2 <- compare_sources(simmr_out_2_combine,
    source_names = c("Z and G", "Enteromorpha"),
    group = 2
  ))
  expect_false(mean(cs_1[[1]]) == mean(cs_2[[1]]))

  # Make sure the plots are different
  cp_1 <- plot(simmr_out_2_combine, type = "histogram", group = 1)
  cp_2 <- plot(simmr_out_2_combine, type = "histogram", group = 2)
  expect_false(cp_1$data$Proportion[1] == cp_2$data$Proportion[1])
})

# Problem with working with tibbles
test_that("simmr combine sources multiple groups", {
  library(readxl)
  path <- system.file("extdata", "geese_data.xls", package = "simmr")
  geese_data <- lapply(excel_sheets(path), read_excel, path = path)

  # Separate the data into parts
  targets <- geese_data[[1]]
  sources <- geese_data[[2]]
  TEFs <- geese_data[[3]]
  concdep <- geese_data[[4]]

  # Load the data into simmr
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

  # Run through simmr
  co(geese_simmr_out <- simmr_mcmc(geese_simmr,
    mcmc_control = list(
      iter = 100,
      burn = 10,
      thin = 1,
      n.chain = 2
    )
  ))

  # Combine sources
  simmr_out_4_combine <- combine_sources(geese_simmr_out,
    to_combine = c(
      "U.lactuca",
      "Enteromorpha"
    ),
    new_source_name = "U.lac and Ent"
  )

  expect_s3_class(simmr_out_4_combine, "simmr_output")
  expect_s3_class(simmr_out_4_combine$input, "simmr_input")
  expect_true(length(simmr_out_4_combine$input$source_names) == 3)
  expect_true(nrow(simmr_out_4_combine$input$correction_means) == 3)
  expect_true(nrow(simmr_out_4_combine$input$source_sds) == 3)
  expect_true(simmr_out_4_combine$input$n_sources == 3)
  expect_true(simmr_out_4_combine$input$n_groups == 8)

  # Check that it still works when you reverse the sources (previously crashed)
  simmr_out_5_combine <- combine_sources(simmr_out_4_combine,
    to_combine = c("Grass", "Zostera"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_5_combine, "simmr_output")
})


# Problem with working with tibbles
test_that("simmr combine sources multiple groups", {
  library(readxl)
  path <- system.file("extdata", "geese_data.xls", package = "simmr")
  geese_data <- lapply(excel_sheets(path), read_excel, path = path)

  # Separate the data into parts
  targets <- geese_data[[1]]
  sources <- geese_data[[2]]
  TEFs <- geese_data[[3]]
  concdep <- geese_data[[4]]



  # Load the data into simmr
  geese_simmr <- simmr_load(
    mixtures = as.matrix(targets[, 1:2]),
    source_names = sources$Sources,
    source_means = as.matrix(sources[, 2:3]),
    source_sds = as.matrix(sources[, 4:5]),
    correction_means = as.matrix(TEFs[, 2:3]),
    correction_sds = as.matrix(TEFs[, 4:5]),
    concentration_means = as.matrix(concdep[, 2:3]),
    group = as.factor(paste("Day", targets$Time))
  )

  # Run through simmr
  co(geese_simmr_out_ffvb <- simmr_ffvb(geese_simmr))

  # Combine sources
  simmr_out_4_combine_ffvb <- combine_sources(geese_simmr_out_ffvb,
    to_combine = c(
      "U.lactuca",
      "Enteromorpha"
    ),
    new_source_name = "U.lac and Ent"
  )

  expect_s3_class(simmr_out_4_combine_ffvb, "simmr_output")
  expect_s3_class(simmr_out_4_combine_ffvb$input, "simmr_input")
  expect_true(length(simmr_out_4_combine_ffvb$input$source_names) == 3)
  expect_true(nrow(simmr_out_4_combine_ffvb$input$correction_means) == 3)
  expect_true(nrow(simmr_out_4_combine_ffvb$input$source_sds) == 3)
  expect_true(simmr_out_4_combine_ffvb$input$n_sources == 3)
  expect_true(simmr_out_4_combine_ffvb$input$n_groups == 8)

  # Check that it still works when you reverse the sources (previously crashed)
  simmr_out_5_combine_ffvb <- combine_sources(simmr_out_4_combine_ffvb,
    to_combine = c("Grass", "Zostera"),
    new_source_name = "Z and G"
  )
  expect_s3_class(simmr_out_5_combine_ffvb, "simmr_output")
})
