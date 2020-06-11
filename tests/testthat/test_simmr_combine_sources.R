context("Test combining of sources") 

library(simmr)

# Single group version first
data("geese_data_day1")

# Taken from the simmr_mcmc example
test_that('simmr combine sources single group', {

  # Load into simmr
  simmr_1 = with(geese_data_day1, 
                 simmr_load(mixtures=mixtures,
                            source_names=source_names,
                            source_means=source_means,
                            source_sds=source_sds,
                            correction_means=correction_means,
                            correction_sds=correction_sds,
                            concentration_means = concentration_means))
  # MCMC run
  simmr_1_out = simmr_mcmc(simmr_1,
                           mcmc_control = list(iter = 100, 
                                               burn = 10, 
                                               thin = 1, 
                                               n.chain = 2))
  
  # Combine two sources
  simmr_out_combine = combine_sources(simmr_1_out,
                                      to_combine=c('Zostera','Grass'),
                                      new_source_name='Z+G')
  expect_s3_class(simmr_out_combine, 'simmr_output')
  expect_s3_class(simmr_out_combine$input, 'simmr_input')
  expect_true(length(simmr_out_combine$input$source_names) == 3)
  expect_true(nrow(simmr_out_combine$input$correction_means) == 3)
  expect_true(nrow(simmr_out_combine$input$source_sds) == 3)
  expect_true(simmr_out_combine$input$n_sources == 3)
  
  expect_output(summary(simmr_out_combine))

})

# Taken from the simmr_mcmc example
test_that('simmr combine sources multiple group', {
  
  data("geese_data")
  
  # Load into simmr
  simmr_2 = with(geese_data, 
                 simmr_load(mixtures=mixtures[,1:2],
                            source_names=source_names,
                            source_means=source_means,
                            source_sds=source_sds,
                            correction_means=correction_means,
                            correction_sds=correction_sds,
                            concentration_means = concentration_means,
                            group = mixtures[,3]))
  # MCMC run
  simmr_2_out = simmr_mcmc(simmr_2,
                           mcmc_control = list(iter = 100, 
                                               burn = 10, 
                                               thin = 1, 
                                               n.chain = 2))
  
  # Combine two sources
  simmr_out_2_combine = combine_sources(simmr_2_out,
                                        to_combine=c('Zostera','Grass'),
                                        new_source_name='Z+G')
  expect_s3_class(simmr_out_2_combine, 'simmr_output')
  expect_s3_class(simmr_out_2_combine$input, 'simmr_input')
  expect_true(length(simmr_out_2_combine$input$source_names) == 3)
  expect_true(nrow(simmr_out_2_combine$input$correction_means) == 3)
  expect_true(nrow(simmr_out_2_combine$input$source_sds) == 3)
  expect_true(simmr_out_2_combine$input$n_sources == 3)
  expect_true(simmr_out_2_combine$input$n_groups == 8)
  
  # Make sure the summaries are different
  summ_1 = summary(simmr_out_2_combine,type='statistics', group = 1)
  summ_2 = summary(simmr_out_2_combine,type='statistics', group = 2)
  expect_false(summ_1$statistics[[1]][1,1] == summ_2$statistics[[1]][1,1])
  
})