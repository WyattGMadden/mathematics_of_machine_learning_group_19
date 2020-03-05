#second attempt to include gibbs sampler

library(ape)
library(phylofactor)

# num= rpois(1,120)
# tree <-rtree(num)
# #determine max groups
# max_groups <- 2*num-2
# 
# #grab clades
# clade1 = getPhyloGroups(tree)[sample(max_groups, 1)][1][[1]][[1]]



#clade2 = getPhyloGroups(tree)[sample(max_groups, 1)][1][[1]][[1]]
#ensure clades aren't monophyletic


if (length(clade1) < 2 ) {next}
#if (length(clade2) < 2 ) {next}

#randomly generate data, then increase size of clade 1 and decrease clade 2                                               

BodySize <-rlnorm(num)
BodySize[clade1] <-rlnorm(length(clade1))*4
#BodySize[clade2] <-rlnorm(length(clade2))/4

#create a data matrix, with the body size, tip labels, and a basis function of all 1s
BodySize = as.data.frame(cbind(BodySize, tree$tip.label, basis = 1, intercept = 1))

#now we generate the missing data.
#these are species we haven't observed data for
num_miss = rpois(1,80)
miss_tip = sample(tree$tip.label, num_miss)

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
  
  product[i] = (0.0001 + sum(is.na(BodySize[grp1,]$BodySize_miss))/length(is.na(BodySize[grp1,]$BodySize)) + sum(is.na(BodySize[grp2,]$BodySize_miss))/length(is.na(BodySize[grp2,]$BodySize)))
  #product[i] =  product[i]
  product[i] =  (unlist(product[i]))^2
  
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
  
  mat <- as.matrix(design.mat[,c('intercept_miss', 'basis_miss')])
  source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/gibbs_sampler.R')
  
    burn_in <- 100
    b0_map_mean <-  beta_samples[(burn_in+1):num_mcmc,1] %>% mean()
    b1_map_mean <-  beta_samples[(burn_in+1):num_mcmc,2] %>% mean()
  
    b0_map_var <- beta_samples[(burn_in+1):num_mcmc,1] %>% var()
    b1_map_var <- beta_samples[(burn_in+1):num_mcmc,2] %>% var()
  
    sigma_mean <- sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% mean()
    sigma_mean <- sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% sd()
    
    t_val_mle <- summary(lm(design.mat[,"body_size_miss"] ~ design.mat[,"basis_miss"]))[4][1]$coefficients[2,3]
    b1_mle <- summary(lm(design.mat[,"body_size_miss"] ~ design.mat[,"basis_miss"]))[4][1]$coefficients[2,1]
    
  #only check if bayes estimate is smaller
  #if (abs(b1_map_mean) > abs(b1_mle) == TRUE) {break}
  
  results[i, ] <- unlist(c(b0_map_mean, 
                         b1_map_mean,
                         b0_map_var,
                         b1_map_var,
                         sigma_mean,
                         sigma_mean,
                         b1_mle,
                         t_val_mle, 
                         (unlist(product[i]))))
  
  #return basis function to normal
  BodySize$basis_miss = BodySize$basis_miss^2
  print(i)
}

which.max(results[,2]/results[,4])
which.max(results[,8])

hist(unlist(product))
#shrunk estimates 
plot(results[,7], results[,2])

#a larger product means more uncertainty in the data, so if the product increases then
#the difference between the MLE and the MAP should increase
# if the product is zero, then the MLE and MAP difference should be zero 
plot(results[,9], abs(results[,7])- abs(results[,2]))

#as the product gets bigger, the shrinking should as well 
lm(abs(results[,7])- abs(results[,2])~ results[,9]) %>% summary()

#need to check if some estimates are larger
which(abs(results[,7])- abs(results[,2]) < 0)

hist(abs(results[,7])- abs(results[,2]))

#plot(results[,1]/results[,2], results[,3]/results[,4])

beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% HPDinterval()
beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% effectiveSize()

tibble(beta0 = beta_samples[(burn_in+1):num_mcmc,1], iteration = 1:(num_mcmc - burn_in)) %>%
ggplot(aes(y = beta0,iteration)) + geom_line()



