#' Compare dietary proportions for a single source across different groups
#' 
#' This function takes in an object of class \code{simmr_output} and creates
#' probabilistic comparisons for a given source and a set of at least two
#' groups.
#' 
#' When two groups are specified, the function produces a direct calculation of
#' the probability that one group is bigger than the other. When more than two
#' groups are given, the function produces a set of most likely probabilistic
#' orderings for each combination of groups. The function produces boxplots by
#' default and also allows for the storage of the output for further analysis if
#' required.
#' 
#' @param simmr_out An object of class \code{simmr_output} created from
#' \code{\link{simmr_mcmc}}.
#' @param source_name The name of a source. This should match the names exactly
#' given to \code{\link{simmr_load}}.
#' @param groups The integer values of the group numbers to be compared. At
#' least two groups must be specified.
#' @param plot A logical value specifying whether plots should be produced or
#' not.
#' 
#' @import ggplot2
#' 
#' @return If there are two groups, a vector containing the differences between
#' the two groups proportions for that source. If there are multiple groups, a
#' list containing the following fields: \item{Ordering }{The different
#' possible orderings of the dietary proportions across groups} \item{out_all
#' }{The dietary proportions for this source across the groups specified as
#' columns in a matrix}
#' @author Andrew Parnell <andrew.parnell@@mu.ie>
#' @seealso See \code{\link{simmr_mcmc}} for complete examples.
#' @examples
#' \dontrun{
#' mix = matrix(c(10.22, 10.37, 10.44, 10.52, 10.19, 10.45, 9.91, 11.27, 
#'                9.34, 11.68, 12.29, 11.04, 11.46, 11.73, 12.29, 11.79, 11.49, 
#'                11.73, 11.1, 11.36, 12.19, 11.03, 11.21, 10.58, 11.61, 12.16, 
#'                10.7, 11.47, 12.07, 11.75, 11.86, 12.33, 12.36, 11.13, 10.92, 
#'                12.42, 10.95, 12.28, 11.04, 10.76, 10.99, 10.78, 11.07, 10.2, 
#'                11.67, 7.53, 10.65, 10.58, 11.13, 7.73, 10.79, 10.47, 10.82, 
#'                10.41, 11.1, 10.95, 10.76, 10.83, 10.25, 10.52, 9.94, 9.94, 11.61, 
#'                10.65, 10.76, 11.11, 10.2, 11.27, 10.21, 10.88, 11.21, 11.36, 
#'                10.75, 12.38, 11.16, 11.57, 10.79, 11.13, 10.72, 10.99, 10.38, 
#'                10.95, 10.75, 10.75, 11.05, 10.66, 10.61, 10.9, 11.14, 10.33, 
#'                10.83, 10.75, 9.18, 9.03, 9.05, 8.6, 8.29, 10.32, 10.28, 6.47, 
#'                11.36, 10.75, 11.13, 11.37, 10.86, 10.54, 10.39, 10.66, 9.99, 
#'                11.65, 11.02, 10.67, 8.15, 11.12, 10.95, 11.2, 10.76, 11.32, 
#'                10.85, 11.74, 10.46, 10.93, 12.3, 10.67, 11.51, 10.56, 12.51, 
#'                13.51, 11.98, 12.2, 10.48, 12.4, 13, 11.36, 12.08, 12.39, 12.28, 
#'                12.6, 11.3, 11.1, 11.42, 11.49, 12, 13.35, 11.97, 13.35, 12.75, 
#'                12.55, 12.3, 12.51, 12.61, 10.98, 11.82, 12.27, 12.11, 12.11, 
#'                12.89, 12.99, 12.29, 11.89, 12.74, 12.29, 11.89, 10.56, 9.27, 
#'                10.54, 10.97, 10.46, 10.56, 10.86, 10.9, 11.06, 10.76, 10.64, 
#'                10.94, 10.85, 10.45, 11.15, 11.23, 11.16, 10.94, 11.2, 10.71, 
#'                9.55, 8.6, 9.67, 8.17, 9.81, 10.94, 9.49, 9.46, 7.94, 9.77, 8.07, 
#'                8.39, 8.95, 9.83, 8.51, 8.86, 7.93, 8, 8.33, 8, 9.39, 8.01, 7.59, 
#'                8.26, 9.49, 8.23, 9.1, 8.21, 9.59, 9.37, 9.47, 8.6, 8.23, 8.39, 
#'                8.24, 8.34, 8.36, 7.22, 7.13, 10.64, 8.06, 8.22, 8.92, 9.35, 
#'                7.32, 7.66, 8.09, 7.3, 7.33, 7.33, 7.36, 7.49, 8.07, 8.84, 7.93, 
#'                7.94, 8.74, 8.26, 9.63, 8.85, 7.55, 10.05, 8.23, 7.74, 9.12, 
#'                7.33, 7.54, 8.8, -11.36, -11.88, -10.6, -11.25, -11.66, -10.41, 
#'                -10.88, -14.73, -11.52, -15.89, -14.79, -17.64, -16.97, -17.25, 
#'                -14.77, -15.67, -15.34, -15.53, -17.27, -15.63, -15.94, -14.88, 
#'                -15.9, -17.11, -14.93, -16.26, -17.5, -16.37, -15.21, -15.43, 
#'                -16.54, -15, -16.41, -15.09, -18.06, -16.27, -15.08, -14.39, 
#'                -21.45, -22.52, -21.25, -21.84, -22.51, -21.97, -20.23, -21.64, 
#'                -22.49, -21.91, -21.65, -21.37, -22.9, -21.13, -19.33, -20.29, 
#'                -20.56, -20.87, -21.07, -21.69, -21.17, -21.74, -22.69, -21.06, 
#'                -20.42, -21.5, -20.15, -21.99, -22.3, -21.71, -22.48, -21.86, 
#'                -21.68, -20.97, -21.91, -19.05, -22.78, -22.36, -22.46, -21.52, 
#'                -21.84, -21.3, -21.39, -22.1, -21.59, -20.14, -20.67, -20.31, 
#'                -20.07, -21.2, -20.44, -22.06, -22.05, -21.44, -21.93, -22.47, 
#'                -22.27, -22.19, -22.81, -20.48, -22.47, -18.06, -20.72, -20.97, 
#'                -19.11, -18.4, -20.45, -21.2, -19.74, -20.48, -21.48, -17.81, 
#'                -19.77, -22.56, -14.72, -12.21, -12.35, -13.88, -14.43, -14.65, 
#'                -13.9, -14.12, -10.88, -10.44, -15.33, -13.78, -13.98, -15.22, 
#'                -15.25, -15.76, -15.78, -15.49, -13.02, -15.3, -15.55, -14.35, 
#'                -14.99, -14.83, -16.18, -15.01, -12.87, -14.67, -13.84, -14.89, 
#'                -13.33, -15.04, -14.29, -15.62, -13.99, -15.06, -15.06, -15, 
#'                -14.55, -13.32, -14.34, -14.47, -14.31, -14.18, -16.18, -16.25, 
#'                -15.92, -15.35, -14.29, -15.92, -15.35, -20.22, -21.4, -19.97, 
#'                -20.78, -20.61, -20.58, -20.19, -20.71, -20.59, -20.09, -19.37, 
#'                -20.41, -20.84, -20.75, -20.29, -20.89, -19.69, -20.41, -21.24, 
#'                -19.33, -25.87, -25.4, -27.23, -27.52, -24.55, -17.36, -24.7, 
#'                -27.76, -28.92, -25.98, -26.77, -28.76, -27.7, -24.75, -25.47, 
#'                -26.58, -28.94, -29.13, -26.65, -28.04, -27.5, -29.28, -27.85, 
#'                -27.41, -27.57, -29.06, -25.98, -28.21, -25.27, -14.43, -27.4, 
#'                -27.76, -28.45, -27.35, -28.83, -29.39, -28.86, -28.61, -29.27, 
#'                -20.32, -28.21, -26.3, -28.27, -27.75, -28.55, -27.38, -29.13, 
#'                -28.66, -29.02, -26.04, -26.06, -28.52, -28.51, -27.93, -29.07, 
#'                -28.41, -26.42, -27.71, -27.75, -24.28, -28.43, -25.94, -28, 
#'                -28.59, -22.61, -27.34, -27.35, -29.14), ncol=2, nrow=251)
#' colnames(mix) = c('d13C','d15N')
#' s_names = c("Zostera", "Grass", "U.lactuca", "Enteromorpha")
#' s_means = matrix(c(6.49, 4.43, 11.19, 9.82, -11.17, -30.88, -11.17, 
#'                    -14.06), ncol=2, nrow=4)
#' s_sds = matrix(c(1.46, 2.27, 1.11, 0.83, 1.21, 0.64, 1.96, 1.17), ncol=2, nrow=4)
#' c_means = matrix(c(3.54, 3.54, 3.54, 3.54, 1.63, 1.63, 1.63, 1.63), ncol=2, nrow=4)
#' c_sds = matrix(c(0.74, 0.74, 0.74, 0.74, 0.63, 0.63, 0.63, 0.63), ncol=2, nrow=4)
#' conc = matrix(c(0.03, 0.04, 0.02, 0.01, 0.36, 0.4, 0.21, 0.18), ncol=2, nrow=4)
#' grp = as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
#'         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
#'         3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 
#'         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
#'         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 
#'         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 
#'         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
#'         7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 
#'         8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8))
#' 
#' # Load this in:
#' simmr_in = simmr_load(mixtures=mix,
#'                      source_names=s_names,
#'                      source_means=s_means,
#'                      source_sds=s_sds,
#'                      correction_means=c_means,
#'                      correction_sds=c_sds,
#'                      concentration_means = conc,
#'                      group=grp)
#' 
#' # Print
#' simmr_in
#' 
#' # Plot
#' plot(simmr_in,group=1:8,xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
#'      ylab=expression(paste(delta^15, "N (\u2030)",sep="")), 
#'      title='Isospace plot of Inger et al Geese data')
#' 
#' # Run MCMC for each group
#' simmr_out = simmr_mcmc(simmr_in)
#' 
#' # Print output
#' simmr_out
#' 
#' # Summarise output
#' summary(simmr_out,type='quantiles',group=1)
#' summary(simmr_out,type='quantiles',group=c(1,3))
#' summary(simmr_out,type=c('quantiles','statistics'),group=c(1,3))
#' 
#' # Plot - only a single group allowed
#' plot(simmr_out,type='boxplot',group=2,title='simmr output group 2')
#' plot(simmr_out,type=c('density','matrix'),grp=6,title='simmr output group 6')
#' 
#' # Compare groups
#' compare_groups(simmr_out,source='Zostera',groups=1:2)
#' compare_groups(simmr_out,source='Zostera',groups=1:3)
#' compare_groups(simmr_out,source='U.lactuca',groups=c(4:5,7,2))
#' }
#' 
#' @export
compare_groups = function(simmr_out,
                          source_name = simmr_out$input$source_names[1],
                          groups = 1:2,
                          plot = TRUE) {
  UseMethod('compare_groups')
}
#' @export
compare_groups.simmr_output = function(simmr_out,
                          source_name = simmr_out$input$source_names[1],
                          groups = 1:2,
                          plot = TRUE) {

# Function to compare between groups both via textual output and with boxplots
# Things to supply are:
# If two groups are given: 
#   - provide the probability of one group being bigger than the other
#   - give the probability distribution of the difference
#   - optional boxplot of two 
# If more than two groups are given:
#   - provide the top most likely orderings of the groups
# An optional boxplot of the groups
  
# Throw an error if only one group is specified
if(length(groups)==1) stop("Please use plot(...) or summary(...) if you just want to look at one group.")
  
# Throw an error if the source name given doesn't match the source names
if(!source_name%in%simmr_out$input$source_names) stop("This source name not found in the current source names. Be sure to check case and spelling")
  
# Get group names
group_names = levels(simmr_out$input$group)[groups]
  
# Start with two groups version
if(length(groups)==2) {
  # Get the output for this particular source on these two groups  
  match_name = match(source_name, simmr_out$input$source_names)
  out_all_grp_1 = simmr_out$output[[groups[1]]]$BUGSoutput$sims.matrix[,match_name]
  out_all_grp_2 = simmr_out$output[[groups[2]]]$BUGSoutput$sims.matrix[,match_name]
  # Produce the difference between the two
  out_diff = out_all_grp_1 - out_all_grp_2

  cat(paste("Prob ( proportion of",source_name,'in group',group_names[1],'> proportion of',source_name,'in group',group_names[2],') =',round(mean(out_diff>0),3)))
  
  if(plot) {
    # Stupid fix for packaging ggplot things
    Group = Proportion = NULL
    df = data.frame(Proportion=c(out_all_grp_1,out_all_grp_2),Group=c(rep(paste(group_names[1]),length(out_all_grp_1)),rep(paste(group_names[2]),length(out_all_grp_2))))
    p = ggplot(df,aes(x=Group,y=Proportion,fill=Group)) + geom_boxplot(alpha=0.5,outlier.size=0) + theme_bw() + theme(legend.position='none') + ggtitle(paste("Comparison of dietary proportions for groups",group_names[1],'and',group_names[2],'for source',source_name))
    print(p)
  }
  
} 

# Now for more groups  
if(length(groups)>2) {
  # Get the output for all the groups
  match_name = match(source_name, simmr_out$input$source_names)
  len = length(simmr_out$output[[groups[1]]]$BUGSoutput$sims.matrix[,match_name])
  out_all = matrix(NA,nrow=len,ncol=length(groups))
  for(j in 1:length(groups)) out_all[,j] = simmr_out$output[[groups[j]]]$BUGSoutput$sims.matrix[,match_name]
  colnames(out_all) = paste(group_names)
  
  # Now find the ordering of each one
  ordering_num = t(apply(out_all,1,order,decreasing=TRUE))
  Ordering = rep(NA,length=nrow(ordering_num))
  for(i in 1:length(Ordering)) Ordering[i] = paste0(group_names[ordering_num[i,]],collapse=" > ")
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
    Group = Proportion = NULL
    df = reshape2::melt(out_all)[,2:3]
    colnames(df) = c('Group','Proportion')
    p = ggplot(df,aes(x=factor(Group),y=Proportion,fill=factor(Group))) + 
      scale_fill_viridis(discrete=TRUE) + 
      geom_boxplot(alpha=0.5,outlier.size=0) + 
      xlab('Group') +
      theme_bw() + theme(legend.position='none') + 
      ggtitle(paste("Comparison of dietary proportions for source",source_name))
    print(p)
  }
  
}  

# Return output
if(length(groups)==2) {
  invisible(list(out_diff))
} else {
  invisible(list(Ordering=Ordering,out_all=out_all))
}  

}
