
library(ape)
library(phylofactor)
set.seed(1)

n= 120
tree <-rtree(n)
plot.phylo(tree,main='Phylogeny - tip labels t1-t5 and tip indexes 1-5')
tiplabels(1:5,cex=1,adj = -2)
edgelabels()

#max groups
getPhyloGroups(tree)[[2*n-2]]

clade1 <- phangorn::Descendants(tree,128,'tips')[[1]]
clade2 <- phangorn::Descendants(tree,186,'tips')[[1]]
clade3 <- phangorn::Descendants(tree,131,'tips')[[1]]
                                                 
set.seed(1)
BodySize <-rlnorm(100)
BodySize[clade1] <-rlnorm(length(clade1))*4
BodySize[clade2] <-rlnorm(length(clade2))/4

BodySize = as.data.frame(cbind(BodySize, tree$tip.label, basis = 1))
miss_tip = sample(tree$tip.label, 20)

#first drop tips from tree
tree_miss <- ape::drop.tip(tree, tree$tip.label[tree$tip.label %in% miss_tip])
BodySize$BodySize_miss = BodySize$BodySize
BodySize$basis_miss = BodySize$basis

BodySize[(BodySize$V2 %in% miss_tip),]$BodySize_miss = NA
BodySize[(BodySize$V2 %in% miss_tip),]$basis_miss = NA
head(BodySize)

BodySize$basis = as.numeric(BodySize$basis)
BodySize$BodySize = as.numeric(BodySize$BodySize)

head(BodySize)

results <- matrix(0, (2*n-2), 5)


product <- list()
for (i in 1: (2*n-2))
     {
  grp1 = getPhyloGroups(tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(tree)[i][[1]][[2]]
  
  #product of two groups, missing observations
  
  product[i] = 0.0001 + sum(is.na(BodySize[grp1,]$BodySize_miss))/length(is.na(BodySize[grp1,]$BodySize)) * sum(is.na(BodySize[grp2,]$BodySize_miss))/length(is.na(BodySize[grp2,]$BodySize))

  #grp1_miss = getPhyloGroups(tree_miss)[i][[1]][[1]]
  #grp2_miss = getPhyloGroups(tree_miss)[i][[1]][[2]]
  
  #BodySize[grp1,]$basis = 1
  BodySize$basis_miss = c(as.numeric(BodySize[grp2,]$basis_miss)*-1, as.numeric(BodySize[grp1,]$basis_miss)*1)
  #BodySize[grp2,]$basis  = as.numeric(BodySize[grp2,]$basis_miss)*-1
  #BodySize[grp2,]$basis_miss  = (as.numeric(BodySize[grp2,]$basis_miss)*-1)
  #BodySize[grp2,]$basis_miss = 0
  
  #Bayesian Posterior Estimate 
  S_o = product[i][[1]]
  m_o = 0
  cov = (BodySize$basis - mean(BodySize$basis)) %*% (BodySize$basis - mean(BodySize$basis))
  
  design.mat = BodySize
  design.mat[is.na(design.mat)] <- 0
  design.mat[(design.mat) == NA] <- 0
  design.mat
  
  S_n = 1/S_o + design.mat$basis_miss %*% design.mat$basis_miss
  m_mle = solve((design.mat$basis_miss) %*% design.mat$basis_miss) %*% t(design.mat$basis_miss) %*% (BodySize$BodySize)
  m_n = 1/S_n * (1/(S_o ) + (design.mat$basis_miss) %*% design.mat$BodySize)
  
  results[i, ] <- unlist(c(m_mle, m_n, S_n, product[i][[1]], i))
}

which.max(results[,2]/results[,3])

plot(results[,1], results[,2])

