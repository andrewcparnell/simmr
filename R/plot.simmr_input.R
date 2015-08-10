plot.simmr_input <-
function(x,tracers=c(1,2),title='Tracers plot',xlab='d13C',ylab='d15N', sigmas=1, group=1, colour=TRUE,...) {

if(length(group)>1) stop("Multiple group plotting not supported.")
  
# Get mixtures to match current group
curr_mix = x$mixtures[x$group==group,,drop=FALSE]
  
# Throw error if plotting only one isotope
if(ncol(curr_mix)==1) stop("This function only works for two or more tracers")

# First get the mean corrected sources and the sd corrected sources
source_means_c = x$source_means + x$correction_means
source_sds_c = x$source_sds^2 + x$correction_sds^2

# Set up data frame for ggplot - have to do it this stupid way because of cran
x2=c(source_means_c[,tracers[1]],curr_mix[,tracers[1]])
y=c(source_means_c[,tracers[2]],curr_mix[,tracers[2]])
x_lower=c(source_means_c[,tracers[1]]-sigmas*source_sds_c[,tracers[1]],curr_mix[,tracers[1]])
x_upper=c(source_means_c[,tracers[1]]+sigmas*source_sds_c[,tracers[1]],curr_mix[,tracers[1]])
y_lower=c(source_means_c[,tracers[2]]-sigmas*source_sds_c[,tracers[2]],curr_mix[,tracers[2]])
y_upper=c(source_means_c[,tracers[2]]+sigmas*source_sds_c[,tracers[2]],curr_mix[,tracers[2]])
Source=c(x$source_names,rep(' Mixtures',nrow(curr_mix)))
size=c(rep(0.8,x$n_sources),rep(0.5,nrow(curr_mix)))
df=data.frame(x=x2,y=y,x_lower,y_lower,x_upper,y_upper,Source,size)

g=ggplot(data=df,aes(x = x,y = y,colour=Source))+theme_bw()+labs(x=xlab,y=ylab,title=title)+geom_errorbarh(aes(xmax=x_upper,xmin=x_lower,height=0))+geom_pointrange(aes(x=x,y=y,ymax=y_upper,ymin=y_lower,height=0.2,shape=Source),size=0.7)+theme(legend.title=element_blank(),legend.key = element_blank())+if(!colour) scale_colour_grey()

print(g)

}
