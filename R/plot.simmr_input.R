plot.simmr_input <-
function(x,tracers=c(1,2),title='Tracers plot',xlab='d13C',ylab='d15N', sigmas=1, colour=TRUE,...) {

# Throw error if plotting only one isotope
if(ncol(x$mixtures)==1) stop("This function only works for two or more tracers")

# First get the mean corrected sources and the sd corrected sources
source_means_c = x$source_means + x$correction_means
source_sds_c = x$source_sds^2 + x$correction_sds^2

# Set up data frame for ggplot - have to do it this stupid way because of cran
x2=c(source_means_c[,tracers[1]],x$mixtures[,tracers[1]])
y=c(source_means_c[,tracers[2]],x$mixtures[,tracers[2]])
x_lower=c(source_means_c[,tracers[1]]-sigmas*source_sds_c[,tracers[1]],x$mixtures[,tracers[1]])
x_upper=c(source_means_c[,tracers[1]]+sigmas*source_sds_c[,tracers[1]],x$mixtures[,tracers[1]])
y_lower=c(source_means_c[,tracers[2]]-sigmas*source_sds_c[,tracers[2]],x$mixtures[,tracers[2]])
y_upper=c(source_means_c[,tracers[2]]+sigmas*source_sds_c[,tracers[2]],x$mixtures[,tracers[2]])
Source=c(x$source_names,rep(' Mixtures',x$n_obs))
size=c(rep(0.8,x$n_sources),rep(0.5,x$n_obs))
df=data.frame(x=x2,y=y,x_lower,y_lower,x_upper,y_upper,Source,size)

g=ggplot(data=df,aes(x = x,y = y,colour=Source))+theme_bw()+labs(x=xlab,y=ylab,title=title)+geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0))+geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source),size=0.7)+theme(legend.title=element_blank(),legend.key = element_blank())+if(!colour) scale_colour_grey()

print(g)

}
