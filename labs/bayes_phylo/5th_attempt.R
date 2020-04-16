
results <- matrix(0, (2*(length(train_tree$tip.label))-2), 9  )

#......................................phylofactor step......................................................
product <- list()
for (i in 1: (2*(length(train_tree$tip.label))-2))
{
  grp1 = getPhyloGroups(train_tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(train_tree)[i][[1]][[2]]
  
  #product of two groups, missing observations
  #as the product gets larger, the estimae goes to zero
  #here a big product means one of the groups is missing a lot of variables
  # ad 0.0001 so we dont get a situation where we divide by zero. bad! 
  product[i] = 1/ (0.0001 + (sum(!is.na(train_BodySize[grp1,]$BodySize_miss))/length(train_BodySize[grp1,]$BodySize) * sum(!is.na(train_BodySize[grp2,]$BodySize_miss))/length(train_BodySize[grp2,]$BodySize)))

  #product[i] =  product[i]
  #product[i] =  (unlist(product[i]))^2
  
  #create constrast basis 
  
  train_BodySize$basis_miss = as.numeric(as.character(train_BodySize$basis_miss))
  train_BodySize$intercept_miss = as.numeric(as.character(train_BodySize$intercept_miss))
  train_BodySize$train_BodySize_miss = as.numeric(as.character(train_BodySize$BodySize_miss))
  
  design.mat = cbind(c(train_BodySize[grp2,]$basis_miss*-1,train_BodySize[grp1,]$basis_miss*1),
                     c(train_BodySize[grp2,]$intercept_miss*1,train_BodySize[grp1,]$intercept_miss*1),
                     c(as.numeric(as.character(train_BodySize[grp2,]$BodySize_miss)),
                       as.numeric(as.character(train_BodySize[grp1,]$BodySize_miss))))
  design.mat[is.na(design.mat)] <- 0
  colnames(design.mat) <- c("basis_miss", "intercept_miss", "body_size_miss")

  #remove the missing data
  design.mat = design.mat[- which(design.mat[,2] == 0),]

  mat <- as.matrix(design.mat[,c('intercept_miss', 'basis_miss')])
  source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/gibbs_sampler.R')
  coef = summary(lm(design.mat[,"body_size_miss"] ~ design.mat[,"basis_miss"]))[4][1]$coefficients
  
  if(dim(coef)[1] > 1)
  {
  t_val_mle <- coef[2,3]
  b1_mle    <- coef[2,1]
  b0_mle    <- coef[1,1]
  }
  if(dim(coef)[1] == 1) # no slope model, no clade split
  {
    t_val_mle <- 0
    b1_mle    <- 0
    b0_mle    <- coef[1,1]
  }
  #only check if bayes estimate is smaller
  #if (abs(b1_map_mean) > abs(b1_mle) == TRUE) {break}
  
  results[i, ] <- unlist(c(b0_map_mean, 
                           b1_map_mean,
                           b0_map_var,
                           b1_map_var,
                           sigma_mean,
                           b0_mle,
                           b1_mle,
                           t_val_mle, 
                           (unlist(product[i]))))
  
  #return basis function to normal
  train_BodySize$basis_miss = train_BodySize$basis_miss^2
  #print(i)
}
#order: 

# b0_map_mean, 
# b1_map_mean,
# b0_map_var,
# b1_map_var,
# sigma_mean,
# b0_mle,
# b1_mle,
# t_val_mle, 
# (unlist(product[i]))
print("done with simulations")

colnames(results) <- c("b0_map_mean", 
                       "b1_map_mean",
                       "b0_map_var",
                       "b1_map_var",
                       "sigma_mean",
                       "b0_mle",
                       "b1_mle",
                       "t_val_mle",
                       "product")
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 

#for the data, select the species and the body size variable 
Data = train_BodySize[,c('Species', 'BodySize_miss')]
#should both be 0:
sum(Data$Species %in% train_tree$tip.label == FALSE)
sum(train_tree$tip.label %in% Data$Species == FALSE)

#convert the variables to the best format
Data <- Data %>% 
  dplyr::mutate(BodySize_miss=  as.numeric(as.character(BodySize_miss))) 

#set the row names of the data to the species, I think this is unecessary now
#remove the NA data, this is the 'missing' data in this simulation
#we haven't observed these tips...sad! 
Data <- Data %>% tidyr::drop_na()
#drop the tips from the tree 

#must covert to a character before dropping, muy importante 
train_tree$tip.label = as.character(train_tree$tip.label)
tree_gpf <- ape::drop.tip(train_tree,train_tree$tip.label[!(train_tree$tip.label %in% Data$Species)])

#okay these should both be 0
sum(Data$Species %in% tree_gpf$tip.label == FALSE)
sum(tree_gpf$tip.label %in% Data$Species == FALSE)

gpf_results = gpf(Data= Data, 
                  tree = tree_gpf, 
                  frmla.phylo = BodySize_miss ~ phylo, algorithm = 'phylo',
                  nfactors = 1)




# rm(beta_samples, coef, cov_beta, design.mat, exp_beta,
#    mat, product, Sigma_0, Sigma_0_inv, X, 
#    b0_map_mean,
#    b0_map_var,
#    b0_mle,
#    b1_map_mean,
#    b1_map_var,
#    b1_mle,
#    beta_0, 
#    burn_in)


# gpf_results$groups
# gpf_results$models
# getPhyloGroups(tree)[samp]
# 
# 
# 
# plot(results[,7], results[,2], ylab = "MLE", xlab = "MAP")
# hist(abs(results[,6]) - abs(results[,1]), xlab = "Intecept MLE-MAP, should be > 0")
# hist(abs(results[,7]) - abs(results[,2]), xlab = "MLE-MAP, should be > 0")
# 
# plot(abs(results[,7]) - abs(results[,2]), (unlist(product)), ylim = c(0,4), xlab = "MLE-MAP, should be > 0")
# 
# hist(abs(results[,8]), xlab = "T Statistic MLE")
# 
# hist(abs(results[,4]/results[,5]), xlab = "MAP Mean / Sigma", breaks = 20)
# 
# plot(results[,4]/results[,5], results[,4], xlab = "MAP Mean / Sigma")
# 
# which.max(abs(results[,8]))
# getPhyloGroups(tree)[which.max(results[,4]/results[,5])][1][[1]][[1]]
# results[samp,8]
# max(results[,8])




