#' Summarises the output created with \code{\link{simmr_mcmc_tdf}}
#' 
#' Produces textual summaries and convergence diagnostics for an object created
#' with \code{\link{simmr_mcmc_tdf}}. The different options are: 'diagnostics'
#' which produces Brooks-Gelman-Rubin diagnostics to assess MCMC convergence,
#' 'quantiles' which produces credible intervals for the parameters,
#' 'statistics' which produces means and standard deviations, and
#' 'correlations' which produces correlations between the parameters.
#' 
#' The quantile output allows easy calculation of 95 per cent credible
#' intervals of the posterior proportions. The Gelman diagnostic values should be close to 1 to ensure satisfactory convergence.
#' 
#' Multiple groups are not currently supported for estimating TDFs.
#' 
#' @param object An object of class \code{simmr_output_df} produced by the
#' function \code{\link{simmr_mcmc_tdf}}
#' @param type The type of output required. At least none of 'diagnostics',
#' 'quantiles', 'statistics', 'correlations'
#' @param ...  Not used
#' @return An list containing the following components: \item{gelman }{The
#' convergence diagnostics} \item{quantiles }{The quantiles of each parameter
#' from the posterior distribution} \item{statistics }{The means and standard
#' deviations of each parameter} \item{correlations }{The posterior
#' correlations between the parameters} Note that this object is reported
#' silently so will be discarded unless the function is called with an object
#' as in the example below.
#' 
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' 
#' @seealso See \code{\link{simmr_mcmc}} for creating objects suitable for this
#' function, and many more examples. See also \code{\link{simmr_load}} for
#' creating simmr objects, \code{\link{plot.simmr_input}} for creating isospace
#' plots, \code{\link{plot.simmr_output}} for plotting output.
#' 
#' @importFrom stats sd cor
#' 
#' @examples
#' \dontrun{
#' ## Example of estimating TDFs for a simple system with known dietary proportions
#' 
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' # Assume p = c(0.25, 0.25, 0.25, 0.25)
#' 
#' # The data
#' mix = matrix(c(-10.13, -10.72, -11.39, -11.18, -10.81, -10.7, -10.54, 
#' -10.48, -9.93, -9.37, 11.59, 11.01, 10.59, 10.97, 11.52, 11.89, 
#' 11.73, 10.89, 11.05, 12.3), ncol=2, nrow=10)
#' colnames(mix) = c('d13C','d15N')
#' s_names=c('Source A','Source B','Source C','Source D')
#' s_means = matrix(c(-14, -15.1, -11.03, -14.44, 3.06, 7.05, 13.72, 5.96), ncol=2, nrow=4)
#' s_sds = matrix(c(0.48, 0.38, 0.48, 0.43, 0.46, 0.39, 0.42, 0.48), ncol=2, nrow=4)
#' conc = matrix(c(0.02, 0.1, 0.12, 0.04, 0.02, 0.1, 0.09, 0.05), ncol=2, nrow=4)
#' 
#' # Load into simmr
#' simmr_tdf = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      concentration_means = conc)
#' 
#' # Plot
#' plot(simmr_tdf)
#' 
#' # MCMC run
#' simmr_tdf_out = simmr_mcmc_tdf(simmr_tdf, 
#' p = matrix(rep(1/simmr_tdf$n_sources, 
#' simmr_tdf$n_sources),
#' ncol = simmr_tdf$n_sources, 
#' nrow = simmr_tdf$n_obs, byrow = TRUE))
#' 
#' # Summary
#' summary(simmr_tdf_out,type='diagnostics')
#' summary(simmr_tdf_out,type='quantiles')
#' 
#' # Now put these corrections back into the model and check the 
#' # iso-space plots and dietary output
#' simmr_tdf_2 = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means = simmr_tdf_out$c_mean_est,
#'                      correction_sds = simmr_tdf_out$c_sd_est,
#'                      concentration_means = conc)
#' 
#' # Plot with corrections now
#' plot(simmr_tdf_2)
#' 
#' simmr_tdf_2_out = simmr_mcmc(simmr_tdf_2)
#' summary(simmr_tdf_2_out, type = 'diagnostics')
#' plot(simmr_tdf_2_out, type = 'boxplot')
#' }
#' 
#' @export
summary.simmr_output_tdf =
  function(object,type=c('diagnostics','quantiles','statistics','correlations'),...) {
    # Get the specified type
    type=match.arg(type,several.ok=TRUE)

    out_all = object$output$BUGSoutput$sims.matrix

    # Get objects
    out_bgr = object$output$BUGSoutput$summary[,'Rhat']
    out_quantiles = t(apply(out_all,2,'quantile',probs=c(0.025,0.25,0.5,0.75,0.975)))
    #  coda:::summary.mcmc.list(object$output)$quantiles
    out_statistics = t(apply(out_all,2,function(x) {return(c(mean=mean(x),sd=stats::sd(x)))}))
    # coda:::summary.mcmc.list(object$output)$statistics[,1:2]
    out_cor = stats::cor(out_all)

    if ('diagnostics'%in%type) {
      # Print out gelman diagnostics of the output
      cat('Gelman diagnostics - these values should all be close to 1.\n')
      cat('If any are larger than 1.1, try a longer run of simmr_mcmc_tdf.\n')
      print(round(out_bgr,2))
    }

    if ('quantiles'%in%type) {
      # Print out quantiles argument
      print(round(out_quantiles,3))
    }

    if ('statistics'%in%type) {
      # Print out quantiles argument
      print(round(out_statistics,3))
    }

    if ('correlations'%in%type) {
      # Print out quantiles argument
      print(round(out_cor,3))
    }

  invisible(list(gelman=out_bgr,quantiles=out_quantiles,statistics=out_statistics,correlations=out_cor))

}
