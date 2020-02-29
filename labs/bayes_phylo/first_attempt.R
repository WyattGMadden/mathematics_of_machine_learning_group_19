
library(ape)
library(phylofactor)
set.seed(1)

#create a tree with n tips
n= 120
tree <-rtree(n)

#plot tree if you would like:
#plot.phylo(tree,main='Phylogeny - tip labels t1-t5 and tip indexes 1-5')
#tiplabels(1:5,cex=1,adj = -2)
#edgelabels()

#determine max groups
max_groups <- 2*n-2

#identify the two clades we would like

clade1 = getPhyloGroups(tree)[104][1][[1]][[1]]
clade2 = getPhyloGroups(tree)[128][1][[1]][[1]]

#clade1 <- phangorn::Descendants(tree,128,'tips')[[1]]
#clade2 <- phangorn::Descendants(tree,186,'tips')[[1]]
#clade3 <- phangorn::Descendants(tree,131,'tips')[[1]]

#randomly generate data, then increase size of clade 1 and decrease clade 2                                               
set.seed(1)
BodySize <-rlnorm(100)
BodySize[clade1] <-rlnorm(length(clade1))*4
BodySize[clade2] <-rlnorm(length(clade2))/4

#create a data matrix, with the body size, tip labels, and a basis function of all 1s
BodySize = as.data.frame(cbind(BodySize, tree$tip.label, basis = 1))

#now we generate the missing data.
#these are species we haven't observed data for
num_miss = 20
miss_tip = sample(tree$tip.label, num_miss)

#first drop tips from tree
#tree_miss <- ape::drop.tip(tree, tree$tip.label[tree$tip.label %in% miss_tip])
BodySize$BodySize_miss = BodySize$BodySize
BodySize$basis_miss = BodySize$basis

#create a new variable, and label the boddy size and basis functions NA for these tips
BodySize[(BodySize$V2 %in% miss_tip),]$BodySize_miss = NA
BodySize[(BodySize$V2 %in% miss_tip),]$basis_miss = NA
head(BodySize)

BodySize$basis = as.numeric(BodySize$basis)
BodySize$BodySize = as.numeric(as.character(BodySize$BodySize))

BodySize$BodySize %>% hist()
head(BodySize)

results <- matrix(0, (2*n-2), 6)

#......................................phylofactor step......................................................
product <- list()
for (i in 1: (2*n-2))
     {
  grp1 = getPhyloGroups(tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(tree)[i][[1]][[2]]
  
  #product of two groups, missing observations
  
  product[i] = (0.0001 + sum(is.na(BodySize[grp1,]$BodySize_miss))/length(is.na(BodySize[grp1,]$BodySize)) * sum(is.na(BodySize[grp2,]$BodySize_miss))/length(is.na(BodySize[grp2,]$BodySize)))

  #create constrast basis 
  BodySize$basis_miss = c(as.numeric(BodySize[grp2,]$basis_miss)*-1, as.numeric(BodySize[grp1,]$basis_miss)*1)
  
  
  #Bayesian Posterior Estimate
  S_o = product[i][[1]]
  m_o = 0

  design.mat = BodySize
  design.mat[is.na(design.mat)] <- 0
  design.mat[(design.mat) == NA] <- 0
  design.mat
  
  #S_n is S_n inverse from the Bishop text, just in case this is confusing to follow
  S_n = 1/S_o + design.mat$basis_miss %*% design.mat$basis_miss #its essentially always 100 plus something extra
  m_mle = solve(design.mat$basis_miss %*% design.mat$basis_miss) %*% t(design.mat$basis_miss) %*% (BodySize$BodySize)
  m_n = 1/S_n * ( design.mat$basis_miss %*% design.mat$BodySize)
  
  if ((abs(m_mle) > abs(m_n)) == FALSE ) {break}
  
  resid <- (sqrt((as.numeric(design.mat$BodySize_miss) - m_mle %*% as.numeric(design.mat$basis_miss))^2))
  var_mle = 1/(n - num_miss) * sum(resid, na.rm = TRUE) %*% solve((design.mat$basis_miss) %*% design.mat$basis_miss)
  
  results[i, ] <- unlist(c(m_mle, var_mle, m_n, S_n, product[i][[1]], i))
  
  #return basis function to normal
  BodySize$basis_miss = BodySize$basis_miss^2
}

which.max(results[,3]/results[,4])
# 64
getPhyloGroups(tree)[64]
which.max(results[,1]/results[,2])
# 64 also
getPhyloGroups(tree)[64]

plot(results[,1]/results[,2], results[,3]/results[,4])

pf_fisher <-twoSampleFactor(log(BodySize$BodySize),tree,nfactors=2,method = "Fisher",ncores = 2)

clade1
clade2
pf_fisher$groups
getPhyloGroups(tree)[2]

ggtree::ggtree(tree) +
   ggtree::geom_hilight(clade1)

phangorn::Descendants(tree,104,'tips')[[1]]

#how to make predictions for new clades 

group = which.max(results[,3]/results[,4])
grp1 = getPhyloGroups(tree)[group][[1]][[1]]
grp2 = getPhyloGroups(tree)[group][[1]][[2]]

#create the predictions based on mu 

#create a new variable, and label the boddy size and basis functions NA for these tips
BodySize$BodySize_test = BodySize$BodySize
BodySize$basis_test = BodySize$basis

BodySize[!(BodySize$V2 %in% miss_tip),]$BodySize_test = NA
BodySize[!(BodySize$V2 %in% miss_tip),]$basis_test = NA


BodySize$basis_test = c(as.numeric(BodySize[grp2,]$basis_test)*-1, as.numeric(BodySize[grp1,]$basis_test)*1)

(BodySize$basis_test * results[group,3] - BodySize$BodySize_test)^2 %>% sum(na.rm=TRUE)
(BodySize$basis_test * results[group,1] - BodySize$BodySize_test)^2 %>% sum(na.rm=TRUE)

