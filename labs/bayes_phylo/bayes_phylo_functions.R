library(ape)
library(phylofactor)
library(tidyverse)


sim_tree <- function(num.tips) {
  num <- num.tips
  tree <-rtree(num)
  tree$tip.label <- as.character(1:num)
  
  #determine max groups
  max_groups <- 2*num-2
  
  #choose a clade of the good length, we want to make sure its not monophyletic 
  clade1 <- 0
  clade2 <- 0
  
  while((length(clade1) < num / 3) | (length(clade2) < num / 3))
  {
    #grab clades
    samp <- sample(max_groups, 1)
    clade1 <- getPhyloGroups(tree)[samp][1][[1]][[1]]
    clade2 <- getPhyloGroups(tree)[samp][1][[1]][[2]]
  }
  
  #randomly generate data                                            
  BodySize <-rlnorm(num, sdlog = 1.5)
  BodySize[clade1] <-rlnorm(length(clade1), sdlog = 1.5)*4
  
  #create a data matrix, with the body size, tip labels, and a basis function of all 1s
  BodySize <- as.data.frame(cbind(BodySize, tree$tip.label, basis = 1, intercept = 1))
  BodySize <- BodySize %>%
    dplyr::mutate(Species=  as.character(V2)) %>%
    dplyr::select(-V2)
  
  #now wegenerate the missing data.
  #these are species we haven't observed data for
  #in this simulation around 1/3 of the data are missing
  
  #num_miss <- 1 + rpois(1,nrow(train_BodySize)/2.5)
  num_miss <- nrow(BodySize)/3
  
  #num_miss <- rpois(1,10)  
  miss_tip <- sample(tree$tip.label, num_miss)
  
  #create missing body size variable
  BodySize$BodySize_miss <- BodySize$BodySize
  BodySize$basis_miss <- BodySize$basis
  BodySize$intercept_miss <- BodySize$intercept
  
  #create a new variable, and label the boddy size and basis functions NA for these tips
  BodySize[(BodySize$Species %in% miss_tip),]$BodySize_miss <- NA
  BodySize[(BodySize$Species %in% miss_tip),]$basis_miss <- NA
  BodySize[(BodySize$Species %in% miss_tip),]$intercept_miss <- NA
  
  BodySize$basis <- as.numeric(BodySize$basis)
  BodySize$BodySize <- as.numeric(as.character(BodySize$BodySize))
  BodySize$intercept <- as.numeric(BodySize$intercept)
  
  tree.data <- BodySize
 
  return(list(tree = tree, 
              tree.data = tree.data, 
              clade.1 = clade1, 
              clade.2 = clade2))
}



drop_tips <- function(tree, tree.data, drop.tips.num) {
  test_tree <- ape::drop.tip(tree,
                             tree$tip.label[!(tree$tip.label %in% sample(tree$tip.label, drop.tips.num))])
  train_tree <- ape::drop.tip(tree,tree$tip.label[(tree$tip.label %in% test_tree$tip.label)])
  
  
  #now, split the dataset into the training and testing datasets
  
  train_BodySize <- tree.data %>% 
    dplyr::filter(Species %in% train_tree$tip.label)
  
  test_BodySize <- tree.data %>% 
    dplyr::filter(Species %in% test_tree$tip.label)
  
  
  train_BodySize <- dplyr::mutate_at(train_BodySize,
                                     c('basis_miss', 'BodySize_miss', 'intercept_miss'),
                                     as.numeric)
  
  test_BodySize <- dplyr::mutate_at(test_BodySize,
                                     c('intercept', 'BodySize', 'basis'),
                                     function(x) as.numeric(as.character(as.numeric(x))))
  
  return(list(tree.test = test_tree,
              tree.train = train_tree, 
              test.data = test_BodySize, 
              train.data = train_BodySize))
}



phylo_bayes <- function(tree.test, tree.train, test.data, train.data, drop.tips.num, epsilon, delta = 1) {
  
  results <- matrix(0, (2*(length(tree.train$tip.label)) - 2), 7)
  colnames(results) <- c("sse_map", "beta_ridge_2", "beta_ridge_1", "N1", "N2", "theta_2", "SSE_train")
  
  #......................................phylofactor step......................................................
  product <- list()
  for (i in 1:(2*(length(tree.train$tip.label)) - 2)) {
    
    grp1 <- getPhyloGroups(tree.train)[i][[1]][[1]]
    grp2 <- getPhyloGroups(tree.train)[i][[1]][[2]]
    
    y <- train.data[,c("BodySize_miss")]
    
    Z <- cbind(c(train.data[grp2,]$basis_miss * - 1, train.data[grp1,]$basis_miss*1),
              c(train.data[grp2,]$intercept_miss * 1, train.data[grp1,]$intercept_miss * 1))
    
    Z <- as.matrix(Z[- which(is.na(Z[,2])),])
    y <- as.numeric(as.character((y[- which(is.na(y))])))
    
    #as theta gets bigger, estimates shrink to 0
    #N1 is proportion  of tips with data, so a higher N1 is better
    #N1 of 0 means every tip in that clade is missing data
    #So as N1 gets bigger, theta 2 should shrink to 0, and our coefficients should go to the LSE
    N1 <- sum(!is.na(train.data[grp1,]$BodySize_miss))/length(train.data[grp1,]$BodySize) 
    N2 <- sum(!is.na(train.data[grp2,]$BodySize_miss))/length(train.data[grp2,]$BodySize)
    
    theta_2 <-  epsilon * exp(-delta * N1 * N2) 
    
    #theta <- diag(2) * c(theta_2, 1/1000)
    theta <- diag(2) * c(theta_2, theta_2)
    
    #theta <- diag(2) * c(0,0)
    beta_ridge <- solve(t(Z) %*% Z + theta) %*% t(Z) %*% y
    #beta <- solve(t(Z) %*% Z) %*% t(Z) %*% y
    
    SSE_train <- 1/nrow(Z) * sum((Z %*% beta_ridge - y)^2)
    #solve(t(Z) %*% Z) %*% t(Z) %*% y
    
    #now cross validate
    grp1_map <- tree.train$tip.label[getPhyloGroups(tree.train)[i][1][[1]][[1]]]
    grp2_map <- tree.train$tip.label[getPhyloGroups(tree.train)[i][1][[1]][[2]]]
    grps <- list(grp1_map, grp2_map)
    cross_val_map <- crossVmap(grps, 
                               tree = tree, 
                               original_community = tree.train$tip.label, 
                               tree.test$tip.label, 
                               ignore.interruptions = F)
    
    cross_val_map <- crossVmap(grps, 
                               tree = tree, 
                               original_community = tree.train$tip.label, 
                               tree.test$tip.label)
    
    #create the basis in the testing set
    basis_map <- c(as.numeric(test.data[test.data$Species %in% cross_val_map[1][[1]],]$basis) * -1, 
                  as.numeric(test.data[test.data$Species %in% cross_val_map[2][[1]],]$basis) * 1)
    
    tryCatch({
      sse_map <- (1 / drop.tips.num ) * sum((((test.data[,"intercept"]*beta_ridge[2] + 
                                                 t(basis_map)*beta_ridge[1]) - test.data[,"BodySize"])^2))}, 
    error = function(e){sse_map <- NA})
    

    
    results[i, 1] <-  sse_map
    results[i, 2] <-  beta_ridge[2]
    results[i, 3] <-  beta_ridge[1]
    results[i, 4] <-  N1
    results[i, 5] <-  N2
    results[i, 6] <-  theta_2
    results[i, 7] <-  SSE_train
    
  }

  
  #for the data, select the species and the body size variable 
  Data <- train.data[,c('Species', 'BodySize_miss')]
  
  #convert the variables to the best format
  Data <- Data %>% 
    dplyr::mutate(BodySize_miss = as.numeric(as.character(BodySize_miss))) 
  
  #set the row names of the data to the species, I think this is unecessary now
  #remove the NA data, this is the 'missing' data in this simulation
  #we haven't observed these tips...sad! 
  Data <- Data %>% tidyr::drop_na()
  #drop the tips from the tree 
  
  #must covert to a character before dropping, muy importante 
  tree.train$tip.label <- as.character(tree.train$tip.label)
  tree_gpf <- ape::drop.tip(tree.train, tree.train$tip.label[!(tree.train$tip.label %in% Data$Species)])
  
  
  gpf_results <- gpf(Data = Data, 
                    tree = tree_gpf, 
                    frmla.phylo = BodySize_miss ~ phylo, algorithm = 'phylo',
                    nfactors = 1)
  
  grp1 <- Data$Species[gpf_results$groups[1][[1]][[1]]]
  grp2 <- Data$Species[gpf_results$groups[1][[1]][[2]]]
  grps <- list(grp1, grp2)
  
  skip_to_next <- FALSE
  
  tryCatch({
    cross_val_gpf <- crossVmap(grps, 
                               tree = tree, 
                               original_community = tree.train$tip.label, 
                               tree.test$tip.label, 
                               ignore.interruptions = F)}, 
    error = function(e){ skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }    
  
  if (length(cross_val_gpf) > 2) {
    basis_grp <-  NA
    sse_gpf <- NA}
  if (length(cross_val_gpf) == 2) {
    basis_grp <- c(as.numeric(test.data[test.data$Species %in% cross_val_gpf[1][[1]],]$basis) * -1, 
                   as.numeric(test.data[test.data$Species %in% cross_val_gpf[2][[1]],]$basis) * 1)
    sse_gpf <- (1 / basis_grp %>% length() ) * sum(((test.data[,"intercept"]*gpf_results$models[[1]]$coefficients[1] + 
                                                       basis_grp*gpf_results$models[[1]]$coefficients[2]) - 
                                                      test.data[,"BodySize"])^2) }
  return(list(results.data = results,
              sse.gpf = sse_gpf))
              
}

phylo_bayes_opti_epsilon <- function(epsilon.seq, iters, tree, tree.data, drop.tips.num, clade.1, clade.2, delta = 1) {
  tries <- 1:iters
  epsilon_seq <- epsilon.seq
  
  error <- matrix(0, length(tries) * 1 * length(epsilon_seq), 22)
  colnames(error) <- c("sse_map_1", "beta_ridge_2_1", "beta_ridge_1_1", "N1_1", "N2_1", "theta_2_1", "SSE_train_1",
                       "sse_map_2", "beta_ridge_2_2", "beta_ridge_1_2", "N1_2", "N2_2", "theta_2_2", "SSE_train_2",
                       "delta", "epsilon", "sse_gpf", "min", "min_2", "yy", "stat1", "stat2")
  
  #plot(epsilon * exp(-delta *seq(0,1,.01)), type = 'l', xlab = "Percent Missing in Tips", ylab ="Penalty")
  #plot(epsilon * exp(-16 *seq(0,1,.01)), type = 'l', xlab = "Percent Missing in Tips", ylab ="Penalty")
  
  
  for(yy in tries) {
    
    dropped_tree_temp <- drop_tips(tree = tree, 
                                   tree.data = tree.data,
                                   drop.tips.num = drop.tips.num)
    tree_train <- dropped_tree_temp[['tree.train']]
    tree_test <- dropped_tree_temp[['tree.test']]
    train_data <- dropped_tree_temp[['train.data']]
    test_data <- dropped_tree_temp[['test.data']]
    
    x <- NULL
    for(x in 1:length(epsilon_seq)) {
      
      epsilon <- epsilon_seq[x]
      
      j <- NULL
      
      #delta is a control parameter for the bayesian model. as it gets bigger the estimates shrink to 0
      #delta = 
      #create the predictions based on mu_map
      #create a new variable, and label the boddy size and basis functions NA for these tip
      #we need to test the max estimate, and the mle estimate
      #when there are no missing data, they find the same split point
      #however, when the tips are missing they tend to get different
      #we should probably test them both on the actual split point
      #first, identify the two groups pulled out by GPF in the training datas
      
      phylo_temp <- phylo_bayes(tree.test = tree_test,
                                tree.train = tree_train,
                                train.data = train_data,
                                test.data = test_data, 
                                drop.tips.num = drop.tips.num,
                                epsilon = epsilon)
      
      results <- phylo_temp[['results.data']]
      sse_gpf <- phylo_temp[['sse.gpf']]
      
      min = which.min(results[, 'sse_map'])
      min_2 = which.min(results[, 'beta_ridge_2'])
      
      
      stat1 = sum(train_data[getPhyloGroups(tree_train)[min][[1]][[1]],]$Species %in% clade.1) / length(clade.1)
      stat2 = sum(train_data[getPhyloGroups(tree_train)[min][[1]][[1]],]$Species %in% clade.2) / length(clade.2)
      
      colnames(results) <- c("sse_map", "beta_ridge_2", "beta_ridge_1", "N1", "N2", "theta_2", "SSE_train")
      iter <- length(epsilon_seq)*(yy - 1) + x
      
      error[iter, 1]   <-  results[min, 'sse_map'] #sse 
      error[iter, 2]   <-  results[min, 'beta_ridge_2']
      error[iter, 3]   <-  results[min, 'beta_ridge_1']
      error[iter, 4]   <-  results[min, 'N1']
      error[iter, 5]   <-  results[min, 'N2']
      error[iter, 6]   <-  results[min, 'theta_2']
      error[iter, 7]   <-  results[min, 'SSE_train']
      
      error[iter, 8]   <-  results[min_2, 'sse_map'] #sse 
      error[iter, 9]   <-  results[min_2, 'beta_ridge_2']
      error[iter, 10]  <-  results[min_2, 'beta_ridge_1']
      error[iter, 11]  <-  results[min_2, 'N1']
      error[iter, 12]  <-  results[min_2, 'N2']
      error[iter, 13]  <-  results[min_2, 'theta_2']
      error[iter, 14]  <-  results[min_2, 'SSE_train']
      
      error[iter, 15]  <-  delta
      error[iter, 16]  <-  epsilon
      error[iter, 17]  <-  sse_gpf
      error[iter, 18]  <-  min
      error[iter, 19]  <-  min_2
      error[iter, 20]  <-  yy
      
      error[iter, 21]  <-  stat1
      error[iter, 22]  <-  stat2
      
      
      print(paste("index = " , iter))
    }
  }
  
  return(error)
}
