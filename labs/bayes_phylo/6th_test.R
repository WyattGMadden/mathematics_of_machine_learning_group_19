library(ape)
library(phylofactor)
sim = (seq(0.01,2,.2)) #must be the same as line 11
eps = (seq(0.01,2,.2))
error <- matrix(0, length(sim)*length(eps), 12)
index = 0
source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/5th_attempt_generate_data.R')

for( x in sim)
{
  print(paste("X =", x))
  epsilon = x
  j=1
  for(j in eps)
  {
    print(paste("j =", j))
    index = index + 1
    delta = 1 #delta is a control parameter for the bayesian model. as it gets bigger the estimates shrink to 0
    delta = delta*(j)
    #delta = delta*j
    #print(num)
    source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/6th_attempt.R')
  
    #create the predictions based on mu_map
    #create a new variable, and label the boddy size and basis functions NA for these tips
  
    #we need to test the max estimate, and the mle estimate
    #when there are no missing data, they find the same split point
    #however, when the tips are missing they tend to get different
    #we should probably test them both on the actual split point..

    #first, identify the two groups pulled out by GPF in the training dataset
  
    #grp1 = train_tree$tip.label[gpf_results$groups[1][[1]][[1]]]
    #grp2 = train_tree$tip.label[gpf_results$groups[1][[1]][[2]]]
    min = which.min(results[,1])
    min_2 = which.min(results[,7])
  
  # results[i,1] <-  sse_map
  # results[i,2] <-  beta_ridge[2]
  # results[i,3] <-  beta_ridge[1]
  # results[i,4] <-  N1
  # results[i,5] <-  N2
  # results[i,6] <-  theta_2
  
   error[index,1]  <-  results[min,1] #sse 
   error[index,2]  <-  results[min,2]
   error[index,3]  <-  results[min,3]
   error[index,4]  <-  results[min,4]
   error[index,5]  <-  results[min,5]
   error[index,6]  <-  results[min,6]
   error[index,7]  <-  delta
   error[index,8]  <-  epsilon
   error[index,9]  <-  sse_gpf
   error[index,10] <-  results[min,7]
   error[index,11] <-  results[min_2,1]
   error[index,12] <-  results[min_2,7]
   print(error[index,])
  }
}

plot((error[,1]), type = 'l', lwd = 5, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.9))
plot((error[,12]), type = 'l', lwd = 5, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.9))   
plot(error[,1] - error[,9], col = rgb(red = 0, green = 0, blue = 1, alpha = 0.9))                         

save(list=ls(),file = "/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/data/test_2020_04_20_A")


