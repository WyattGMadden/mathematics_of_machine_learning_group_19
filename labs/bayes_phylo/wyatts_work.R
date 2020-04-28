library(phylofactor)
library(ape)

set.seed(Sys.time())
num= 40 #rpois(1,60)
tree <-rtree(num)
tree$tip.label <- as.character(1:num)
drop_tips = num/2

#determine max groups
max_groups <- 2*num-2

#choose a clade of the good length, we want to make sure its not monophyletic 
clade1 = 0
clade2 = 0

while((length(clade1) < num / 3) | (length(clade2) < num / 3)) {
  #grab clades
  samp = sample(max_groups, 1)
  clade1 = getPhyloGroups(tree)[samp][1][[1]][[1]]
  clade2 = getPhyloGroups(tree)[samp][1][[1]][[2]]
}

#randomly generate data                                            
BodySize <-rlnorm(num)
BodySize[clade1] <-rlnorm(length(clade1))*6

#create a data matrix, with the body size, tip labels, and a basis function of all 1s
BodySize = as.data.frame(cbind(BodySize, tree$tip.label, basis = 1, intercept = 1))
BodySize <- BodySize %>%
  dplyr::mutate(Species=  as.character(V2)) %>%
  dplyr::select(-V2)

#grab a sample of 10 tips for the testing tree
test_tree <- ape::drop.tip(tree, sample(tree$tip.label, drop_tips))
train_tree <- ape::drop.tip(tree, test_tree$tip.label)

#should be equal to 0
sum(test_tree$tip.label %in% train_tree$tip.label) == 0
sum(train_tree$tip.label %in% test_tree$tip.label) == 0
num - (length(train_tree$tip.label) + length(test_tree$tip.label)) == 0

#now, split the dataset into the training and testing datasets

train_BodySize = BodySize %>% 
  dplyr::filter(Species %in% train_tree$tip.label)

test_BodySize = BodySize %>% 
  dplyr::filter(Species %in% test_tree$tip.label)

#should both be true:
sum(train_BodySize$Species %in% test_BodySize$Species) == 0
sum(test_BodySize$Species %in% train_BodySize$Species) == 0

#now wegenerate the missing data.
#these are species we haven't observed data for
#in this simulation around 1/3 of the data are missing

#num_miss = 1 + rpois(1,nrow(train_BodySize)/2.5)
num_miss = nrow(train_BodySize)/2
#num_miss = rpois(1,10)  
miss_tip = sample(train_tree$tip.label, num_miss)

#create missing body size variable
train_BodySize$BodySize_miss = train_BodySize$BodySize
train_BodySize$basis_miss = train_BodySize$basis
train_BodySize$intercept_miss = train_BodySize$intercept

#create a new variable, and label the boddy size and basis functions NA for these tips
train_BodySize[(train_BodySize$Species %in% miss_tip),]$BodySize_miss = NA
train_BodySize[(train_BodySize$Species %in% miss_tip),]$basis_miss = NA
train_BodySize[(train_BodySize$Species %in% miss_tip),]$intercept_miss = NA

train_BodySize$basis = as.numeric(train_BodySize$basis)
train_BodySize$BodySize = as.numeric(as.character(train_BodySize$BodySize))
train_BodySize$intercept = as.numeric(train_BodySize$intercept)





test_BodySize <- test_BodySize %>%
  dplyr::mutate(intercept = as.numeric(as.character(intercept))) %>%
  dplyr::mutate(BodySize = as.numeric(as.character(BodySize))) %>%
  dplyr::mutate(basis = as.numeric(as.character(basis)))

train_BodySize$basis_miss = as.numeric(as.character(train_BodySize$basis_miss))
train_BodySize$intercept_miss = as.numeric(as.character(train_BodySize$intercept_miss))
train_BodySize$train_BodySize_miss = as.numeric(as.character(train_BodySize$BodySize_miss))

results <- matrix(0, (2*(length(train_tree$tip.label))-2), 7 )

#......................................phylofactor step......................................................
product <- list()
for (i in 1: (2*(length(train_tree$tip.label))-2)){
  #print(i)
  grp1 = getPhyloGroups(train_tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(train_tree)[i][[1]][[2]]
  
  y = train_BodySize[,c("BodySize_miss")]
  
  Z = cbind(c(train_BodySize[grp2,]$basis_miss*-1,train_BodySize[grp1,]$basis_miss*1),
            c(train_BodySize[grp2,]$intercept_miss*1,train_BodySize[grp1,]$intercept_miss*1))
  
  Z = as.matrix(Z[- which(is.na(Z[,2])),])
  y = as.numeric(as.character((y[- which(is.na(y))])))
  
  #as theta gets bigger, estimates shrink to 0
  
  N1 = sum(!is.na(train_BodySize[grp1,]$BodySize_miss))/length(train_BodySize[grp1,]$BodySize) 
  N2 = sum(!is.na(train_BodySize[grp2,]$BodySize_miss))/length(train_BodySize[grp2,]$BodySize)
  
  theta_2 =  epsilon * exp(delta * N1 * N2) 
  
  theta = diag(2) * c(theta_2, 1/1000)
  #theta = diag(2) * c(0,0)
  theta
  beta_ridge = solve(t(Z) %*% Z + theta) %*% t(Z) %*% y
  
  cbind(Z %*% beta_ridge, y)
  
  SSE_train = 1/nrow(Z) * sum((Z %*% beta_ridge - y)^2)
  #solve(t(Z) %*% Z) %*% t(Z) %*% y
  
  #now cross validate
  grp1_map = train_tree$tip.label[getPhyloGroups(train_tree)[i][1][[1]][[1]]]
  grp2_map = train_tree$tip.label[getPhyloGroups(train_tree)[i][1][[1]][[2]]]
  grps=list(grp1_map, grp2_map)
  cross_val_map = crossVmap(grps, tree= tree, original_community = train_tree$tip.label, test_tree$tip.label, ignore.interruptions = F)
  
  cross_val_map = crossVmap(grps, tree= tree, original_community = train_tree$tip.label, test_tree$tip.label)
  
  #if (length(cross_val_map) > 2) next 
  
  #print(length(cross_val_map[1][[1]]))
  #print(length(cross_val_map[2][[1]]))
  
  #create the basis in the testing set
  basis_map = c(as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_map[1][[1]],]$basis) * -1, as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_map[2][[1]],]$basis) * 1)
  
  tryCatch({
    sse_map = (1 / drop_tips ) * sum((((test_BodySize[,"intercept"]*beta_ridge[2] + t(basis_map)*beta_ridge[1]) - test_BodySize[,"BodySize"])^2))
  }, 
  error=function(e){ sse_map <<- NA})
  
  #print(sse_map)
  
  results[i,1] <-  sse_map
  results[i,2] <-  beta_ridge[2]
  results[i,3] <-  beta_ridge[1]
  results[i,4] <-  N1
  results[i,5] <-  N2
  results[i,6] <-  theta_2
  results[i,7] <-  SSE_train
  
}

plot(results[,1])

plot(results[,1], results[,7])

which.min(results[,1])
which.min(results[,7])

getPhyloGroups(train_tree)[which.min(results[,1])]
getPhyloGroups(train_tree)[which.min(results[,7])]
clade1
clade2
#min(results[,1], na.rm=TRUE)


#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
#now run gpr 
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



grp1 = Data$Species[gpf_results$groups[1][[1]][[1]]]
grp2 = Data$Species[gpf_results$groups[1][[1]][[2]]]
grps = list(grp1, grp2)

skip_to_next = FALSE
tryCatch({
  cross_val_gpf = crossVmap(grps, tree= tree, original_community = train_tree$tip.label, test_tree$tip.label, ignore.interruptions = F)
}, error=function(e){ skip_to_next <<- TRUE})

if(skip_to_next) { next }    

if (length(cross_val_gpf) > 2)  
{basis_grp =  NA
sse_gpf = NA
}
if (length(cross_val_gpf) == 2)  
{
  basis_grp = c(as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_gpf[1][[1]],]$basis) * -1, as.numeric(test_BodySize[test_BodySize$Species %in% cross_val_gpf[2][[1]],]$basis) * 1)
  sse_gpf = (1 / basis_grp %>% length() ) * sum(((test_BodySize[,"intercept"]*gpf_results$models[[1]]$coefficients[1] + basis_grp*gpf_results$models[[1]]$coefficients[2]) - test_BodySize[,"BodySize"])^2)
}










