#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................


set.seed(10262018)
y = as.numeric(design.mat[,"body_size_miss"])
X=mat

burn_in = 500

#X <- model.matrix(winpercent ~ chocolate_factor + nut_factor, data = candy)
p <- ncol(X)
n <- sum(X[,1])

# Initialization and Prior
num_mcmc <- 5000
beta_0 <- rep(0,p)

#in Sigma_0, as the covariance gets larger, the estimate approaches the MLE 
#we would like the opposite, as the product gets larger, we would like the estimate to approach 0
#so, lets take the inverse of the product here
#but let the covariance on the intercept be large, let that float 

#product of two groups, missing observations
#as the product gets larger, the estimae goes to zero
#here a big product means one of the groups is missing a lot of variables
# ad 0.0001 so we dont get a situation where we divide by zero. bad! 

N1 = sum(!is.na(train_BodySize[grp1,]$BodySize_miss))/length(train_BodySize[grp1,]$BodySize) 
N2 = sum(!is.na(train_BodySize[grp2,]$BodySize_miss))/length(train_BodySize[grp2,]$BodySize)

#as alpha gets larga and larga, W should shrink to 0 
#if everything is present we want alpha to be 0, least squares 
#plot(y = epsilon * exp(-delta(0,1,.1)), x=  seq(0,1,.1), type = 'l')
#two tuning parameters, delta, and epsilon

#epsilon = 
product[i] = exp(-delta * N1 * N2) 
alpha =  epsilon * exp(-delta * N1 * N2) 
Sigma_0 <- diag(p) * c(100, 1/alpha)
#Sigma_0 <- diag(p) * c(1000, 1/1000)

#we want: product is large, MAP approaches 0 
#small product, map approaches the MLE 

Sigma_0_inv <- solve(Sigma_0)
nu_0 <- .02
sigmasq_0 <- 1
beta_samples <- matrix(0, nrow = num_mcmc, ncol = p)
sigmasq_samples <- rep(1, num_mcmc)

for (iter in 2:num_mcmc){
  # sample beta
  cov_beta <- solve(Sigma_0_inv + t(X) %*% X / sigmasq_samples[iter - 1])
  #exp_beta <- cov_beta %*% (Sigma_0_inv %*% beta_0 + t(X) %*% y / sigmasq_samples[iter-1])
  exp_beta <- cov_beta %*% (t(X) %*% y / sigmasq_samples[iter-1])
  
  beta_samples[iter,] <- mnormt::rmnorm(1, exp_beta, cov_beta) 
  
  # sample sigmasq
  sigmasq_samples[iter] <- LearnBayes::rigamma(1, .5 * (nu_0 + n) , 
                                   .5 * (nu_0 * sigmasq_0 + t(y - X %*% beta_samples[iter,]) %*% 
                                           (y - X %*% beta_samples[iter,])) )
}

b0_map_mean <-  beta_samples[(burn_in+1):num_mcmc,1] %>% mean()
b1_map_mean <-  beta_samples[(burn_in+1):num_mcmc,2] %>% mean()

b0_map_var <- beta_samples[(burn_in+1):num_mcmc,1] %>% var()
b1_map_var <- beta_samples[(burn_in+1):num_mcmc,2] %>% var()

sigma_mean <- sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% mean()
sigma_mean <- sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% sd()

#beta_samples[(burn_in+1):num_mcmc,2] %>% mean()
#beta_samples[(burn_in+1):num_mcmc,2] %>% hist()
#beta_samples[(burn_in+1):num_mcmc,1] %>% hist()

#beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% HPDinterval()
#beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% effectiveSize()

#tibble(beta0 = beta_samples[(burn_in+1):num_mcmc,1], iteration = 1:(num_mcmc - burn_in)) %>% 
#  ggplot(aes(y = beta0,iteration)) + geom_line()

