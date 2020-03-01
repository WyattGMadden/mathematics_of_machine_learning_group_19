sim = 400
error <- matrix(0, sim, 7)

for(j in 1:sim)
{
  source('/Users/dancrowley/Documents/machine_learning_zosso/mathematics_of_machine_learning_group_19/labs/bayes_phylo/first_attempt.R')
  #error group 1

  
  error[j,1] <-  map_error
  #error group 2
  error[j,2] <-  mle_error
  error[j,3] <-  chance_error
  
  error[j,4] <-  n
  error[j,5] <-  num_miss
  error[j,6] <-  group_map
  error[j,7] <-  group_mle
  
  print(j)
  rm(list=setdiff(ls(), c("error", "j")))
}

#bayesian MAP
plot(error[,1], type = 'l')
#red is the MLE
lines(error[,2], type = 'l', col = 'red')
lines(error[,3], type = 'l', col = 'blue')

hist(error[,1] - error[,2], breaks = 20)
hist(error[,1] - error[,3], breaks = 20)
hist(error[,2] - error[,3], breaks = 20)

sum(error[,1] - error[,2] >0)
sum(error[,1] - error[,2] <0)

plot(error[,1], error[,2])


