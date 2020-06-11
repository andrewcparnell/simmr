#' Compare dietary proportions between multiple sources
#' 
#' This function takes in an object of class \code{simmr_output} and creates
#' probabilistic comparisons between the supplied sources. The group number can
#' also be specified.
#' 
#' When two sources are specified, the function produces a direct calculation
#' of the probability that the dietary proportion for one source is bigger than
#' the other. When more than two sources are given, the function produces a set
#' of most likely probabilistic orderings for each combination of sources. The
#' function produces boxplots by default and also allows for the storage of the
#' output for further analysis if required.
#' 
#' @param simmr_out An object of class \code{simmr_output} created from
#' \code{\link{simmr_mcmc}}.
#' @param source_names The names of at least two sources. These should match
#' the names exactly given to \code{\link{simmr_load}}.
#' @param group The integer values of the group numbers to be compared. If not
#' specified assumes the first or only group
#' @param plot A logical value specifying whether plots should be produced or
#' not.
#' 
#' @import ggplot2
#' @importFrom reshape2 "melt"
#' 
#' @return If there are two sources, a vector containing the differences
#' between the two dietary proportion proportions for these two sources. If
#' there are multiple sources, a list containing the following fields:
#' \item{Ordering }{The different possible orderings of the dietary proportions
#' across sources} \item{out_all }{The dietary proportions for these sources
#' specified as columns in a matrix}
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} for complete examples.
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
#' plot(simmr_1_out,type='boxplot')
#' plot(simmr_1_out,type='histogram')
#' plot(simmr_1_out,type='density')
#' plot(simmr_1_out,type='matrix')
#' 
#' # Compare two sources
#' compare_sources(simmr_1_out,source_names=c('Source B','Source D'))
#' 
#' # Compare multiple sources
#' compare_sources(simmr_1_out)
#' }
#' 
#' @export
compare_sources = function(simmr_out,
                           source_names = simmr_out$input$source_names,
                           group = 1,
                           plot = TRUE) {
  UseMethod('compare_sources')
}  
#' @export
compare_sources.simmr_output = function(simmr_out,
                           source_names = simmr_out$input$source_names,
                           group = 1,
                           plot = TRUE) {

# Function to compare between sources within a group both via textual output and with boxplots
# Things to ly are:
# If two sources are given: 
#   - provide the probability of one group having higher dietary proportion than the other
#   - give the probability distribution of the difference
#   - optional boxplot of two 
# If more than two sources are given:
#   - provide the top most likely orderings of the sources
# An optional boxplot of the sources
  
# Throw an error if only one group is specified
if(length(source_names)==1) stop("Use compare_between_groups if you want to compare a single source between groups.")
  
# Throw an error if the source name given doesn't match the source names
if(!all(source_names%in%simmr_out$input$source_names)) stop("Some source names not found in the current source names. Be sure to check case and spelling")
  
# Start with two groups version
if(length(source_names)==2) {
  # Get the output for this particular source on these two groups  
  match_names = match(source_names, simmr_out$input$source_names)
  out_all_src_1 = simmr_out$output[[group]]$BUGSoutput$sims.list$p[,match_names[1]]
  out_all_src_2 = simmr_out$output[[group]]$BUGSoutput$sims.list$p[,match_names[2]]
  # Produce the difference between the two
  out_diff = out_all_src_1 - out_all_src_2
  cat(paste("Prob ( proportion of",source_names[1],'> proportion of',source_names[2],') =',round(mean(out_diff>0),3)))
  
  if(plot) {
    # Stupid fix for packaging ggplot things
    Source = Proportion = NULL
    df = data.frame(Proportion=c(out_all_src_1,out_all_src_2),Source=c(rep(source_names[1],length(out_all_src_1)),rep(source_names[2],length(out_all_src_2))))
    p = ggplot(df,aes(x=Source,y=Proportion,fill=Source)) + geom_boxplot(alpha=0.5,outlier.size=0) + theme_bw() + theme(legend.position='none') + ggtitle(paste("Comparison of dietary proportions for sources",source_names[1],'and',source_names[2]))
    print(p)
  }
  
} 

# Now for more groups  
if(length(source_names)>2) {
  # Get the output for all the groups
  match_names = match(source_names, simmr_out$input$source_names)
  out_all = simmr_out$output[[group]]$BUGSoutput$sims.list$p[,match_names]
  
  # Now find the ordering of each one
  ordering_num = t(apply(out_all,1,order,decreasing=TRUE))
  Ordering = rep(NA,length=nrow(ordering_num))
  for(i in 1:length(Ordering)) Ordering[i] = paste0(source_names[ordering_num[i,]],collapse=" > ")
  cat('Most popular orderings are as follows:\n')
  tab = t(t(sort(table(Ordering,dnn=NULL),decreasing=TRUE)))
  colnames(tab) = 'Probability'
  # Do not print all of it if too long
  if(nrow(tab)>30) {
    print(round(tab[1:30,]/length(Ordering),4))
  } else {
    print(round(tab/length(Ordering),4))
  }
    
  if(plot) {
    # Stupid fix for packaging ggplot things
    Source = Proportion = NULL
    df = reshape2::melt(out_all)[,2:3]
    colnames(df) = c('Source','Proportion')
    p = ggplot(df,aes(x=Source,y=Proportion,fill=Source)) + 
      scale_fill_viridis(discrete=TRUE) + 
      geom_boxplot(alpha=0.5,outlier.size=0) + 
      theme_bw() + 
      theme(legend.position='none') + 
      ggtitle(paste("Comparison of dietary proportions between sources"))
    print(p)
  }
  
}  

# Return output
if(length(source_names)==2) {
  invisible(list(out_diff))
} else {
  invisible(list(Ordering=Ordering,out_all=out_all))
}  

}
