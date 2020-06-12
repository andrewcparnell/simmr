context('Plot input data') 

data('geese_data_day1')
simmr_1 = with(geese_data_day1, 
                simmr_load(mixtures=mixtures,
                           source_names=source_names,
                           source_means=source_means,
                           source_sds=source_sds,
                           correction_means=correction_means,
                           correction_sds=correction_sds,
                           concentration_means = concentration_means))

test_that("basic simmr_input plot", {
  expect_s3_class(plot(simmr_1), 'ggplot')
  expect_s3_class(plot(simmr_1, colour = FALSE), 'ggplot')
  expect_s3_class(plot(simmr_1,tracers=c(2,1)), 'ggplot')
})

test_that('1D simmr plot', {
  simmr_2 = with(geese_data_day1, 
                 simmr_load(mixtures=mixtures[1,,drop=FALSE],
                            source_names=source_names,
                            source_means=source_means,
                            source_sds=source_sds,
                            correction_means=correction_means,
                            correction_sds=correction_sds,
                            concentration_means = concentration_means))
  
  
  # Plot 3 times - first default d13C vs d15N 
  expect_s3_class(plot(simmr_2), 'ggplot')
  expect_s3_class(plot(simmr_2, colour = FALSE), 'ggplot')
})

test_that('Multi-groups plot', {
  data(geese_data)
  simmr_4 = with(geese_data, 
                      simmr_load(mixtures=mixtures,
                                 source_names=source_names,
                                 source_means=source_means,
                                 source_sds=source_sds,
                                 correction_means=correction_means,
                                 correction_sds=correction_sds,
                                 concentration_means = concentration_means,
                                 group=groups))

  expect_s3_class(plot(simmr_4), 'ggplot')
  expect_s3_class(plot(simmr_4, group = 5), 'ggplot')
})

test_that('Single iso plot', {
  data(geese_data)
  simmr_5 = with(geese_data, 
                 simmr_load(mixtures=mixtures[,1,drop = FALSE],
                            source_names=source_names,
                            source_means=source_means[,1,drop = FALSE],
                            source_sds=source_sds[,1,drop = FALSE],
                            correction_means=correction_means[,1,drop = FALSE],
                            correction_sds=correction_sds[,1,drop = FALSE],
                            concentration_means = concentration_means[,1,drop = FALSE],
                            group=groups))
  
  expect_s3_class(plot(simmr_5), 'ggplot')
  expect_s3_class(plot(simmr_5, group = 5), 'ggplot')
})


