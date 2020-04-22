# delta = 1
# epsilon = 1

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
for (i in 1: (2*(length(train_tree$tip.label))-2))
{
  #print(i)
  grp1 = getPhyloGroups(train_tree)[i][[1]][[1]]
  grp2 = getPhyloGroups(train_tree)[i][[1]][[2]]
  
  y = train_BodySize[,c("BodySize_miss")]
  
  Z = cbind(c(train_BodySize[grp2,]$basis_miss*-1,train_BodySize[grp1,]$basis_miss*1),
        c(train_BodySize[grp2,]$intercept_miss*1,train_BodySize[grp1,]$intercept_miss*1))
  
  Z = as.matrix(Z[- which(is.na(Z[,2])),])
  y = as.numeric(as.character((y[- which(is.na(y))])))
  
  #as theta gets bigger, estimates shrink to 0
  #N1 is proportion  of tips with data, so a higher N1 is better
  #N1 of 0 means every tip in that clade is missing data
  #So as N1 gets bigger, theta 2 should shrink to 0, and our coefficients should go to the LSE
  N1 = sum(!is.na(train_BodySize[grp1,]$BodySize_miss))/length(train_BodySize[grp1,]$BodySize) 
  N2 = sum(!is.na(train_BodySize[grp2,]$BodySize_miss))/length(train_BodySize[grp2,]$BodySize)
  
  theta_2 =  epsilon * exp(-delta * N1 * N2) 
  
  #theta = diag(2) * c(theta_2, 1/1000)
  theta = diag(2) * c(theta_2, theta_2)
  
  #theta = diag(2) * c(0,0)
  theta
  beta_ridge = solve(t(Z) %*% Z + theta) %*% t(Z) %*% y
  #beta = solve(t(Z) %*% Z) %*% t(Z) %*% y
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

gpf_results$groups

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









