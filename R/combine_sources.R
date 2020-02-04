#' Combine the proportions from two food sources after running simmr
#'
#' This function takes in an object of class \code{simmr_output} and combines
#' two of the food sources. It works for single and multiple group data.
#'
#' Often two sources either (1) lie in similar location on the iso-space plot,
#' or (2) are very similar in phylogenetic terms. In case (1) it is common to
#' experience high (negative) posterior correlations between the sources.
#' Combining them can reduce this correlation and improve precision of the
#' estimates. In case (2) we might wish to determine the joint amount eaten of
#' the two sources when combined. This function thus combines two sources after
#' a run of \code{\link{simmr_mcmc}} (known as a posteriori combination). The
#' new object can then be called with \code{\link{plot.simmr_input}} or
#' \code{\link{plot.simmr_output}} to produce iso-space plots of summaries of
#' the output after combination.
#'
#' @param simmr_out An object of class \code{simmr_output} created from
#' \code{\link{simmr_mcmc}}.
#' @param to_combine The names of exactly two sources. These should match the
#' names given to \code{\link{simmr_load}}.
#' @param new_source_name A name to give to the new combined source
#' @return A new \code{simmr_output} object
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} and the associated vignette for
#' examples.
#' @examples
#' \dontrun{
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#'
#' # The data
#' mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54,
#' -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89,
#' 11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10)
#' colnames(mix) = c('d13C','d15N')
#' s_names=c('Source A','Source B','Source C','Source D')
#' s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, nrow=4)
#' s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4)
#' c_means = matrix(c(2.63, 1.59, 3.41, 3.04, 3.28, 2.34, 2.14, 2.36), ncol=2, nrow=4)
#' c_sds = matrix(c(0.41, 0.44, 0.34, 0.46, 0.46, 0.48, 0.46, 0.66), ncol=2, nrow=4)
#' conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4)
#'
#' # Load into simmr
#' simmr_1 = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means=c_means,
#'                      correction_sds=c_sds,
#'                      concentration_means = conc)
#'
#' # Plot
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # MCMC run
#' simmr_1_out = simmr_mcmc(simmr_1)
#'
#' # Print it
#' print(simmr_1_out)
#'
#' # Summary
#' summary(simmr_1_out)
#' summary(simmr_1_out,type='diagnostics')
#' summary(simmr_1_out,type='correlations')
#' summary(simmr_1_out,type='statistics')
#' ans = summary(simmr_1_out,type=c('quantiles','statistics'))
#'
#' # Plot
#' plot(simmr_1_out)
#' plot(simmr_1_out,type='boxplot')
#' plot(simmr_1_out,type='histogram')
#' plot(simmr_1_out,type='density')
#' plot(simmr_1_out,type='matrix')
#'
#' simmr_out_combine = combine_sources(simmr_1_out,
#' to_combine=c('Source A','Source D'),
#' new_source_name='Source A+D')
#' plot(simmr_out_combine$input)
#' plot(simmr_out_combine,type='boxplot',title='simmr output: combined sources')
#'
#' }
#'
#' @export 
combine_sources = function(simmr_out,
                           to_combine = simmr_out$input$source_names[1:2],
                           new_source_name = 'combined_source') {
   UseMethod('combine_sources')
}
#' @export 
combine_sources.simmr_output = function(simmr_out,
                           to_combine = simmr_out$input$source_names[1:2],
                           new_source_name = 'combined_source') {
  # A posteriori combining of sources
  
  # Check only two sources to be combined
  if (length(to_combine) != 2)
    stop("Currently only two sources can be combined")
  
  # # Check class
  # if (class(simmr_out) != 'simmr_output')
  #   stop("Only objects of class simmr_output can be run through this function")
  
  # Find which columns to combine by number
  to_combine_cols = match(to_combine, simmr_out$input$source_names)
  if (any(is.na(to_combine_cols)))
    stop('1 or more source names not found')
  
  simmr_new_out = simmr_out
  
  # 1 combine the chosen source means
  old_source_means = simmr_out$input$source_means
  simmr_new_out$input$source_means = old_source_means[-to_combine_cols[2], ,
                                                      drop = FALSE]
  simmr_new_out$input$source_means[to_combine_cols[1], ] = apply(old_source_means[to_combine_cols, ,drop = FALSE], 2, 'mean')
  
  # 2 combine the source sds
  old_source_sds = simmr_out$input$source_sds
  simmr_new_out$input$source_sds = old_source_sds[-to_combine_cols[2], ,drop= FALSE]
  simmr_new_out$input$source_sds[to_combine_cols[1], ] = apply(old_source_sds[to_combine_cols, , drop = FALSE], 2, function(x)
    sqrt(sum(x ^ 2)))
  
  # 3 combine the correction means
  old_correction_means = simmr_out$input$correction_means
  simmr_new_out$input$correction_means = old_correction_means[-to_combine_cols[2], ,drop = FALSE]
  simmr_new_out$input$correction_means[to_combine_cols[1], ] = apply(old_correction_means[to_combine_cols, ,drop = FALSE], 2, 'mean')
  
  # 4 combine the correction sds
  old_correction_sds = simmr_out$input$correction_sds
  simmr_new_out$input$correction_sds = old_correction_sds[-to_combine_cols[2], ,drop = FALSE]
  simmr_new_out$input$correction_sds[to_combine_cols[1], ] = apply(old_correction_sds[to_combine_cols, , drop = FALSE], 2, function(x)
    sqrt(sum(x ^ 2)))
  
  # 5 combine the concentraion means
  old_concentration_means = simmr_out$input$concentration_means
  simmr_new_out$input$concentration_means = old_concentration_means[-to_combine_cols[2], ,drop = FALSE]
  simmr_new_out$input$concentration_means[to_combine_cols[1], ] = apply(old_concentration_means[to_combine_cols, ,drop = FALSE], 2, 'mean')
  
  # 6 change the source names
  old_source_names = simmr_out$input$source_names
  simmr_new_out$input$source_names = old_source_names[-to_combine_cols[2]]
  simmr_new_out$input$source_names[to_combine_cols[1]] = new_source_name
  
  # 7 Change n_sources
  simmr_new_out$input$n_sources = simmr_new_out$input$n_sources - 1
  
  # 8 Sum across all the output values
  n_groups = simmr_out$input$n_groups
  for (j in 1:n_groups) {
    simmr_new_out$output[[j]] = simmr_out$output[[j]]
    # Change sims.list and sims.matrix
    # First sims.matrix
    sims.matrix = simmr_out$output[[j]]$BUGSoutput$sims.matrix
    new_sims.matrix = sims.matrix[,-(to_combine_cols[2]+1)]
    new_sims.matrix[,to_combine_cols[1]+1] = sims.matrix[,to_combine_cols[1]+1] + sims.matrix[,to_combine_cols[2]+1]
    colnames(new_sims.matrix)[to_combine_cols[1]+1] = new_source_name
    simmr_new_out$output[[j]]$BUGSoutput$sims.matrix = new_sims.matrix
    # Now sims.list
    sims.list = simmr_out$output[[1]]$BUGSoutput$sims.list
    new_sims.list = sims.list
    new_sims.list$p = sims.list$p[,-to_combine_cols[2]]
    new_sims.list$p[,to_combine_cols[1]] = sims.list$p[,to_combine_cols[2]] + sims.list$p[,to_combine_cols[1]]
    simmr_new_out$output[[j]]$BUGSoutput$sims.list = new_sims.list
  }
  
  class(simmr_new_out) = 'simmr_output'
  return(simmr_new_out)
  
}
