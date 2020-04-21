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

while((length(clade1) < num / 3) | (length(clade2) < num / 3))
{
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
test_tree <- ape::drop.tip(tree,tree$tip.label[!(tree$tip.label %in% sample(tree$tip.label, drop_tips))])
train_tree <- ape::drop.tip(tree,tree$tip.label[(tree$tip.label %in% test_tree$tip.label)])

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
