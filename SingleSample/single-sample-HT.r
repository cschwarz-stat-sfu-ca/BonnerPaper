# Do a comparison of the expectation and variance of a HT approach

# N   = size of population
# p_i = X_i sampled from beta(alpha,beta) for each unit
# Y_i ~ dbern(pi)

# HT approach
#  condition on Y_i=1
#  generate U*_i = geometric(p_i_
#  estimate N* =- sum(U*) + n


library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)


cat(file="model.txt", "
    ############################################################
model {
    
# input data are n (number of samples), x_i=p_i

  # generate the HT equivalent estimator
  for(i in 1:n){
     m.star[i] ~ dnegbin(x[i],1)
  }
  U.star= sum(m.star)
  N.star= U.star + n
}
    ") # End of the model


# initialize the seed
set.seed(23433)

# generate the data
N <- 10000
alpha = 2
beta  = 2

sim.data <- data.frame(X = rbeta(N,alpha, beta))
sim.data$Y <- rbinom(N, 1, sim.data$X)

head(sim.data)
mean(sim.data$X)
sum(sim.data$Y)

# Next create the data.txt file.

# The datalist will be passed to JAGS with the names of the data
# values.
data.list <- list(n=sum(sim.data$Y),
                  x=sim.data$X[sim.data$Y==1])


# check the list
data.list


init.list <- list(
  list(),
  list(),
  list()
)  # end of list of lists of initial values



# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("U.star","N.star") # parameters to monitor


results <- jags( 
  data      =data.list,   # list of data variables
  inits     =init.list,   # list/function for initial values
  parameters=monitor.list,# list of parameters to monitor
  model.file="model.txt",  # file with bugs model
  n.chains=3,
  n.iter  =5000,          # total iterations INCLUDING burn in
  n.burnin=2000,          # number of burning iterations
  n.thin=2,               # how much to thin
  DIC=FALSE,               # is DIC to be computed?
  working.dir=getwd()    # store results in current working directory
)


# now results is a BIG list of stuff
names(results)
names(results$BUGSoutput)


# get the summary table
results$BUGSoutput$summary
results$BUGSoutput$summary[,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]



