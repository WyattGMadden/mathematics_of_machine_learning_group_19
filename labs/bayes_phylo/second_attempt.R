#second attempt to include gibbs sampler

library(ape)
library(phylofactor)

num= rpois(1,120)
tree <-rtree(num)
#determine max groups
max_groups <- 2*num-2

#grab clades
clade1 = getPhyloGroups(tree)[sample(max_groups, 1)][1][[1]][[1]]
clade2 = getPhyloGroups(tree)[sample(max_groups, 1)][1][[1]][[1]]
#ensure clades aren't monophyletic
if (length(clade1) < 2 ) {next}
if (length(clade2) < 2 ) {next}

#randomly generate data, then increase size of clade 1 and decrease clade 2                                               

BodySize <-rlnorm(num)
BodySize[clade1] <-rlnorm(length(clade1))*4
BodySize[clade2] <-rlnorm(length(clade2))/4

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

results <- matrix(0, (2*num-2), 8  )

#......................................phylofactor step......................................................
product <- list()
for (i in 1: (2*num-2))
{
  grp1 = getPhyloGroups(tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(tree)[i][[1]][[2]]
  
  #product of two groups, missing observations
  
  product[i] = (0.0001 + sum(is.na(BodySize[grp1,]$BodySize_miss))/length(is.na(BodySize[grp1,]$BodySize)) * sum(is.na(BodySize[grp2,]$BodySize_miss))/length(is.na(BodySize[grp2,]$BodySize)))
  #product[i] =  product[i]
  product[i] =  sqrt(unlist(product[i]))
  
  #create constrast basis 
  BodySize$basis_miss = c(as.numeric(BodySize[grp2,]$basis_miss)*-1, as.numeric(BodySize[grp1,]$basis_miss)*1)
  BodySize$intercept_miss = as.numeric(BodySize$intercept_miss)
  
  #Bayesian Posterior Estimate
  S_o = product[i][[1]]
  m_o = 0
  
  design.mat = BodySize
  design.mat$BodySize_miss = as.numeric(design.mat$BodySize_miss)
  design.mat[is.na(design.mat)] <- 0
  
  mat <- as.matrix(design.mat[,c('intercept_miss', 'basis_miss')])
  source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/gibbs_sampler.R')
  
    burn_in <- 100
    b0_map_mean <-  beta_samples[(burn_in+1):num_mcmc,1] %>% mean()
    b1_map_mean <-  beta_samples[(burn_in+1):num_mcmc,2] %>% mean()
  
    b0_map_var <- beta_samples[(burn_in+1):num_mcmc,1] %>% var()
    b1_map_var <- beta_samples[(burn_in+1):num_mcmc,2] %>% var()
  
    sigma_mean <- sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% mean()
    sigma_mean <- sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% sd()
    
    t_val_mle <- summary(lm(design.mat$BodySize_miss ~ design.mat$basis_miss))[4][1]$coefficients[2,3]
    b1_mle <- summary(lm(design.mat$BodySize_miss ~ design.mat$basis_miss))[4][1]$coefficients[2,1]
    
  #only check if bayes estimate is smaller
  #if (abs(b1_map_mean) > abs(b1_mle) == TRUE) {break}
  
  results[i, ] <- unlist(c(b0_map_mean, 
                         b1_map_mean,
                         b0_map_var,
                         b1_map_var,
                         sigma_mean,
                         sigma_mean,
                         b1_mle,
                         t_val_mle))
  
  #return basis function to normal
  BodySize$basis_miss = BodySize$basis_miss^2
  print(i)
}

which.max(results[,2]/results[,4])
which.max(results[,8])


hist(results[,2])

#shrunk estimates 
plot(results[,8], results[,2])

getPhyloGroups(tree)[64]
which.max(results[,1]/results[,2])
getPhyloGroups(tree)[64]

#plot(results[,1]/results[,2], results[,3]/results[,4])

#pf_fisher <-twoSampleFactor(log(BodySize$BodySize),tree,nfactors=2,method = "Fisher",ncores = 2)
#pf_fisher$groups
#getPhyloGroups(tree)[2]

#ggtree::ggtree(tree) +
#   ggtree::geom_hilight(clade1)

#phangorn::Descendants(tree,104,'tips')[[1]]

#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................
#..........................................cross validating................................................

#how to make predictions for new clades 

group_map = which.max(results[,3]/results[,4])
grp1 = getPhyloGroups(tree)[group_map][[1]][[1]]
grp2 = getPhyloGroups(tree)[group_map][[1]][[2]]

#create the predictions based on mu_map
#create a new variable, and label the boddy size and basis functions NA for these tips
BodySize$BodySize_test = BodySize$BodySize
BodySize$basis_test = BodySize$basis
BodySize[!(BodySize$V2 %in% miss_tip),]$BodySize_test = NA
BodySize[!(BodySize$V2 %in% miss_tip),]$basis_test = NA
BodySize$basis_test = c(as.numeric(BodySize[grp2,]$basis_test)*-1, as.numeric(BodySize[grp1,]$basis_test)*1)
#error group 1
map_error = 1/n*(BodySize$basis_test * results[group_map,3] - BodySize$BodySize_test)^2 %>% sum(na.rm=TRUE)
#error group 2

group_mle = which.max(results[,1]/results[,2])
grp1 = getPhyloGroups(tree)[group_mle][[1]][[1]]
grp2 = getPhyloGroups(tree)[group_mle][[1]][[2]]

#create the predictions based on mu_map
#create a new variable, and label the boddy size and basis functions NA for these tips
BodySize$BodySize_test = BodySize$BodySize
BodySize$basis_test = BodySize$basis
BodySize[!(BodySize$V2 %in% miss_tip),]$BodySize_test = NA
BodySize[!(BodySize$V2 %in% miss_tip),]$basis_test = NA
BodySize$basis_test = c(as.numeric(BodySize[grp2,]$basis_test)*-1, as.numeric(BodySize[grp1,]$basis_test)*1)

mle_error = 1/n*(BodySize$basis_test * results[group_mle,1] - BodySize$BodySize_test)^2 %>% sum(na.rm=TRUE)

chance_error = 1/n*(mean(BodySize$BodySize_test, na.rm =TRUE) - BodySize$BodySize_test)^2 %>% sum(na.rm=TRUE)
