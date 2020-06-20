context('simmr_load')

library(simmr)

test_that('simplest example', {
  
  data("geese_data_day1")

  # Load in with simmr_load
  simmr_1 = with(geese_data_day1, 
                 simmr_load(mixtures=mixtures,
                            source_names=source_names,
                            source_means=source_means,
                            source_sds=source_sds,
                            correction_means=correction_means,
                            correction_sds=correction_sds,
                            concentration_means = concentration_means))
  expect_s3_class(simmr_1, 'simmr_input')
  expect_true(is.matrix(simmr_1$source_means))
  expect_true(is.matrix(simmr_1$source_sds))
  expect_true(is.matrix(simmr_1$correction_means))
  expect_true(is.matrix(simmr_1$correction_sds))
  expect_true(is.matrix(simmr_1$concentration_means))

})

test_that('group example', {
  data("geese_data")
  simmr_groups = with(geese_data, 
                      simmr_load(mixtures=mixtures,
                                 source_names=source_names,
                                 source_means=source_means,
                                 source_sds=source_sds,
                                 correction_means=correction_means,
                                 correction_sds=correction_sds,
                                 concentration_means = concentration_means,
                                 group=groups))

  expect_s3_class(simmr_groups, 'simmr_input')
  expect_true(is.matrix(simmr_groups$source_means))
  expect_true(is.matrix(simmr_groups$source_sds))
  expect_true(is.matrix(simmr_groups$correction_means))
  expect_true(is.matrix(simmr_groups$correction_sds))
  expect_true(is.matrix(simmr_groups$concentration_means))
  expect_true(is.factor(simmr_groups$group))
  expect_true(is.integer(simmr_groups$group_int))
})

# Would be good to test for a load of errors at this point

test_that('test without tefs and concentration', {
  
  data("geese_data_day1")
  
  # Load in with simmr_load
  simmr_2 = with(geese_data_day1, 
                 simmr_load(mixtures=mixtures,
                            source_names=source_names,
                            source_means=source_means,
                            source_sds=source_sds))
  expect_s3_class(simmr_2, 'simmr_input')
  expect_true(is.matrix(simmr_2$source_means))
  expect_true(is.matrix(simmr_2$source_sds))
  expect_true(is.matrix(simmr_2$correction_means))
  expect_true(is.matrix(simmr_2$correction_sds))
  expect_true(is.matrix(simmr_2$concentration_means))
  
})