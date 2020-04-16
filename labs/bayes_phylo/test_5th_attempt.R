library(ape)
library(phylofactor)

sim = 100
error <- matrix(0, sim, 12)

for(j in 1:sim)
{
  source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/5th_attempt_generate_data.R')
  print(num)
  source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/5th_attempt.R')
  
  #create the predictions based on mu_map
  #create a new variable, and label the boddy size and basis functions NA for these tips
  
  #we need to test the max estimate, and the mle estimate
  #when there are no missing data, they find the same split point
  #however, when the tips are missing they tend to get different
  #we should probably test them both on the actual split point..
  
  
  #first, identify the two groups pulled out by GPF in the training dataset
 
  grp1 = train_tree$tip.label[gpf_results$groups[1][[1]][[1]]]
  grp2 = train_tree$tip.label[gpf_results$groups[1][[1]][[2]]]
  grps = list(grp1, grp2)
  
  plot.phylo(test_tree,use.edge.length=FALSE);edgelabels()
  
  cross_val_gpf = crossVmap(grps, tree= tree, original_community = train_tree$tip.label, test_tree$tip.label, ignore.interruptions = F)

  if (length(cross_val_gpf) > 2) next 
  

  test_BodySize$basis_grp = c(as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_gpf[1][[1]],]$basis) * -1, as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_gpf[2][[1]],]$basis) * 1)
  
  max_map = which.max(results[,4]/results[,5])
  getPhyloGroups(train_tree)[max_map]
  grp1_map = train_tree$tip.label[getPhyloGroups(train_tree)[max_map][1][[1]][[1]]]
  grp2_map = train_tree$tip.label[getPhyloGroups(train_tree)[max_map][1][[1]][[2]]]
  grps_map = list(grp1_map, grp2_map)
  cross_val_map = crossVmap(grps, tree= tree, original_community = train_tree$tip.label, test_tree$tip.label, ignore.interruptions = F)
  test_BodySize$basis_map = c(as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_map[1][[1]],]$basis) * -1, as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_map[2][[1]],]$basis) * 1)
  
  test_BodySize <- test_BodySize %>%
    dplyr::mutate(intercept = as.numeric(as.character(intercept))) %>%
    dplyr::mutate(BodySize = as.numeric(as.character(BodySize))) %>%
    dplyr::mutate(basis = as.numeric(as.character(basis)))
    
  

  sse_map = (1 / 10 ) * sum(((test_BodySize[,"intercept"]*results[max_map,"b0_map_mean"]          + t(test_BodySize[,"basis_map"])*results[max_map,"b1_map_mean"])          - test_BodySize[,"BodySize"])^2)
  sse_gpf = (1 / 10 ) * sum(((test_BodySize[,"intercept"]*gpf_results$models[[1]]$coefficients[1] + t(test_BodySize[,"basis_grp"])*gpf_results$models[[1]]$coefficients[2]) - test_BodySize[,"BodySize"])^2)
 
   N = sum(!is.na(train_BodySize[,"intercept_miss"]))
  sse_map_train = 1/N*sum(((train_BodySize[,"intercept_miss"]*results[max_map,"b0_map_mean"]          + t(train_BodySize[,"basis_miss"])*results[max_map,"b1_map_mean"])            - as.numeric(as.character(train_BodySize[,"BodySize_miss"])))^2, na.rm=TRUE)
  sse_gpf_train = 1/N*sum(((train_BodySize[,"intercept_miss"]*gpf_results$models[[1]]$coefficients[1] + t(train_BodySize[,"basis_miss"])*gpf_results$models[[1]]$coefficients[2])    - as.numeric(as.character(train_BodySize[,"BodySize_miss"])))^2, na.rm=TRUE)
  
  error[j,1]  <-  sse_gpf
  error[j,2]  <-  sse_map
  error[j,3]  <-  sse_gpf_train
  error[j,4]  <-  sse_map_train
  error[j,5]  <-  num_miss
  error[j,6]  <-  num
  error[j,7]  <-  #max_mle
  error[j,8]  <-  samp
  error[j,9]  <-  NA #results[samp,"b0_mle"]
  error[j,10] <-  NA #results[samp,"b1_mle"]
  error[j,11] <-  NA #results[samp,"b0_map_mean"]
  error[j,12] <-  NA #results[samp,"b1_map_mean"]

  rm(list=setdiff(ls(), c("error", "sim", "j")))
  
  print(j)
}
  
lines((error[,1]), type = 'l', lwd = 5)  #mle testing
plot((error[,2]), col = 'red', lwd = 5, type = 'l')       #map testing
lines((error[,3]), col = 'blue', lwd = 5)                  #mle training
lines((error[,4]), col = 'green', lwd = 5)                 #map training

hist(error[,1] - error[,2], breaks = 10)
hist(error[,4] - error[,3], breaks = 10)

save(list=ls(),file = "/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/data/test_2020_04_15")


