simmr_mcmc = function(simmr_in, prior.control=list(means=rep(0,simmr_in$n_sources),sd=rep(1,simmr_in$n_sources)), mcmc.control=list(iter=10000,burn=1000,thin=10,n.chain=4)) {

# Main function to run simmr through JAGS
if(class(simmr_in)!='simmr_input') stop("Input argument simmr_in must have come from simmr_load")

# Throw warning if n.chain =1
if(mcmc.control$n.chain==1) warning("Running only 1 MCMC chain will cause an error in the convergence diagnostics")

# Determine if a single observation or not
if(nrow(simmr_in$mixtures)==1) {
  cat('Only 1 mixture value, performing a simmr solo run...\n')
  solo=TRUE
} else {
  solo=FALSE
}
  
# Set up the model string
model_string = '
model {
  # Likelihood
  for (j in 1:J) {
    for (i in 1:N) {  
      y[i,j] ~ dnorm(inprod(p*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p,q[,j]), 1/var_y[j])
    }
    var_y[j] <- inprod(pow(p*q[,j],2),pow(s_sd[,j],2)+pow(c_sd[,j],2))/pow(inprod(p,q[,j]),2)
+ pow(sigma[j],2)
  }

  # Prior on sigma
  for(j in 1:J) { sigma[j] ~ dunif(0,sig_upp) }

  # CLR prior on p
  p[1:K] <- expf/sum(expf)
  for(k in 1:K) {
    expf[k] <- exp(f[k])
    f[k] ~ dnorm(mu_f[k],1/pow(sigma_f[k],2))
  }
}
'

# Create data object
data = with(simmr_in,list(
  y=mixtures,
  s_mean=source_means,
  s_sd=source_sds,
  N=n_obs,
  J=n_tracers,
  c_mean=correction_means,
  c_sd = correction_sds,
  q=concentration_means,
  K=n_sources,
  mu_f=prior.control$means,
  sigma_f=prior.control$sd,
  sig_upp=ifelse(solo,0.001,1000)))

# Run in JAGS
model = rjags::jags.model(textConnection(model_string), data=data, n.chain=mcmc.control$n.chain, n.adapt=mcmc.control$burn)
output = rjags::coda.samples(model=model, variable.names=c("p","sigma"), n.iter=mcmc.control$iter, thin=mcmc.control$thin)

output_2 = lapply(output,"colnames<-",c(simmr_in$source_names, paste0('sd_',colnames(simmr_in$mixtures))))
class(output_2) = c('mcmc.list')

output_all = vector('list')
output_all$input = simmr_in
output_all$output = output_2
class(output_all) = 'simmr_output'

return(output_all)

}
