#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................
#....................................gibbs sampler.........................................

#bridger first
#to create our posterior distribution we need to replicate our in class R script 

#### 1. Gibbs Sampler

mat

set.seed(10262018)
y = as.numeric(design.mat$BodySize_miss)
X=mat

#X <- model.matrix(winpercent ~ chocolate_factor + nut_factor, data = candy)
n = ncol(X)
p <- ncol(X)
n <- sum(X[,1])

# Initialization and Prior
num_mcmc <- 5000
beta_0 <- rep(0,p)
Sigma_0 <- diag(p) * product[i][[1]]
Sigma_0_inv <- solve(Sigma_0)
nu_0 <- .02
sigmasq_0 <- 1
beta_samples <- matrix(0, nrow = num_mcmc, ncol = p)
sigmasq_samples <- rep(1, num_mcmc)

for (iter in 2:num_mcmc){
  # sample beta
  cov_beta <- solve(Sigma_0_inv + t(X) %*% X / sigmasq_samples[iter - 1])
  exp_beta <- cov_beta %*% (Sigma_0_inv %*% beta_0 + t(X) %*% y / sigmasq_samples[iter-1])
  beta_samples[iter,] <- rmnorm(1, exp_beta, cov_beta) 
  
  # sample sigmasq
  sigmasq_samples[iter] <- rigamma(1, .5 * (nu_0 + n) , 
                                   .5 * (nu_0 * sigmasq_0 + t(y - X %*% beta_samples[iter,]) %*% 
                                           (y - X %*% beta_samples[iter,])) )
}

burn_in <- 100
beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% summary()
beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% HPDinterval()
beta_samples[(burn_in+1):num_mcmc,] %>% as.mcmc() %>% effectiveSize()

sqrt(sigmasq_samples[(burn_in+1):num_mcmc]) %>% as.mcmc() %>% summary()

tibble(beta0 = beta_samples[(burn_in+1):num_mcmc,1], iteration = 1:(num_mcmc - burn_in)) %>% 
  ggplot(aes(y = beta0,iteration)) + geom_line()

