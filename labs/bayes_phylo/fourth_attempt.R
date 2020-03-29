

#randomly generate data, then increase size of clade 1 and decrease clade 2                                               
BodySize <-rlnorm(num)
BodySize[clade1] <-rlnorm(length(clade1))*4

#create a data matrix, with the body size, tip labels, and a basis function of all 1s
BodySize = as.data.frame(cbind(BodySize, tree$tip.label, basis = 1, intercept = 1))

#first drop tips from tree
#tree_miss <- ape::drop.tip(tree, tree$tip.label[tree$tip.label %in% miss_tip])
BodySize$BodySize_miss = BodySize$BodySize
BodySize$basis_miss = BodySize$basis
BodySize$intercept_miss = BodySize$intercept

#create a new variable, and label the boddy size and basis functions NA for these tips
BodySize[(BodySize$V2 %in% miss_tip),]$BodySize_miss = NA
BodySize[(BodySize$V2 %in% miss_tip),]$basis_miss = NA
BodySize[(BodySize$V2 %in% miss_tip),]$intercept_miss = NA

BodySize$basis = as.numeric(BodySize$basis)
BodySize$BodySize = as.numeric(as.character(BodySize$BodySize))
BodySize$intercept = as.numeric(BodySize$intercept)

results <- matrix(0, (2*num-2), 9  )

#......................................phylofactor step......................................................
product <- list()
for (i in 1: (2*num-2))
{
  grp1 = getPhyloGroups(tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(tree)[i][[1]][[2]]
  
  #product of two groups, missing observations
  #as the product gets larger, the estimae goes to zero
  #here a big product means one of the groups is missing a lot of variables
  # ad 0.0001 so we dont get a situation where we divide by zero. bad! 
  product[i] = 1/ (0.0001 + (sum(!is.na(BodySize[grp1,]$BodySize_miss))/length(BodySize[grp1,]$BodySize) * sum(!is.na(BodySize[grp2,]$BodySize_miss))/length(BodySize[grp2,]$BodySize)))

  #product[i] =  product[i]
  #product[i] =  (unlist(product[i]))^2
  
  #create constrast basis 
  
  BodySize$basis_miss = as.numeric(as.character(BodySize$basis_miss))
  BodySize$intercept_miss = as.numeric(as.character(BodySize$intercept_miss))
  BodySize$BodySize_miss = as.numeric(as.character(BodySize$BodySize_miss))
  
  design.mat = cbind(c(BodySize[grp2,]$basis_miss*-1,BodySize[grp1,]$basis_miss*1),
                     c(BodySize[grp2,]$intercept_miss*1,BodySize[grp1,]$intercept_miss*1),
                     c(BodySize[grp2,]$BodySize_miss,BodySize[grp1,]$BodySize_miss))
  design.mat[is.na(design.mat)] <- 0
  colnames(design.mat) <- c("basis_miss", "intercept_miss", "body_size_miss")
  
  # BodySize$basis_miss = c(as.numeric(as.character(BodySize[grp2,]$basis_miss))*-1, as.numeric(as.character(BodySize[grp1,]$basis_miss))*1)
  # BodySize$intercept_miss = as.numeric(BodySize$intercept_miss)
  
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
  BodySize$basis_miss = BodySize$basis_miss^2
  print(i)
}

plot(results[,7], results[,2], ylab = "MLE", xlab = "MAP")
hist(abs(results[,6]) - abs(results[,1]), xlab = "Intecept MLE-MAP, should be > 0")
hist(abs(results[,7]) - abs(results[,2]), xlab = "MLE-MAP, should be > 0")

plot(abs(results[,7]) - abs(results[,2]), (unlist(product)), ylim = c(0,4), xlab = "MLE-MAP, should be > 0")

hist(abs(results[,8]), xlab = "T Statistic MLE")

which.max(abs(results[,8]))
getPhyloGroups(tree)[which.max(results[,8])][1][[1]][[1]]
results[samp,8]
max(results[,8])

