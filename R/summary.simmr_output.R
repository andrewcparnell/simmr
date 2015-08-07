summary.simmr_output <-
function(object,type=c('diagnostics','quantiles','statistics','correlations'),...) {

  # Get the specified type
  type=match.arg(type,several.ok=TRUE)

  out_all = do.call(rbind,object$output)
  
  out_bgr = coda::gelman.diag(object$output,multivariate=FALSE)$psrf
  out_quantiles = t(apply(out_all,2,'quantile',probs=c(0.025,0.25,0.5,0.75,0.975)))
  #  coda:::summary.mcmc.list(object$output)$quantiles
  out_statistics = t(apply(out_all,2,function(x) {return(c(mean=mean(x),sd=stats::sd(x)))}))
  # coda:::summary.mcmc.list(object$output)$statistics[,1:2]
  out_cor = stats::cor(out_all)
  
  if ('diagnostics'%in%type) {
    # Print out gelman diagnostics of the output
    cat('Gelman diagnostics - these values should all be close to 1.\n')
    cat('If not, try a longer run of simmr_mcmc.\n')
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
