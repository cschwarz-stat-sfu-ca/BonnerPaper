# Implement the simple binomial model using the R2jags template
#
# 2017-02-02 CJS Update
# 2015-02-27 CJS Typo corrections and updates
# 2015-02-09 CJS First edition
#
# This script creates the model.txt file, the data.txt, and inits.txt
# file. It then calls JAGS.
#
# Ensure that the R2jags package is available.

library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)
library(tidyverse)

# The BUGS model is specified as a text file.

# The model file.

cat(file="model.txt", "
############################################################

model {
# Prior distributions for model parameters
    psi   ~ dbeta(.0001,1)
    mu.wt ~ dnorm(0,.001)
    tau.wt~ dgamma(.001,.001)
    sigma <- sqrt(1/tau.wt)
    a0 ~ dnorm(-1,10)
    a1 ~ dnorm( 1,10)

    for(i in 1:(nind+nz)){
       wt[i] ~ dnorm(mu.wt,tau.wt)
       z[i] ~  dbin(psi,1)
       logit(p[i]) <- a0 + a1*wt[i]
       p.obs[i] <- (1-(1-p[i])^2)*z[i]
       obs[i] ~dbin(p.obs[i],1)
    }
    N<-sum(z[1:(nind+nz)])
}
") # End of the model



# Next create the data.txt file.

N <- 1000 # Population size
mu_x <- 0                               # Covariate mean
sigma_x  <- 1                           # Covaraite sd

beta <- c(-1,1)                          # Parameters of capture probability
## Set seed
set.seed(8888)

sim_data <- tibble(x = rnorm(N,mu_x,sigma_x),
                   p = (1 + exp(-(beta[1] + beta[2]*x)))^-1,
                   y1 = rbinom(N,1,p),
                   y2 = rbinom(N,1,p),
                   w = y1 + y2,
                   p.obs= 1-(1-p)^2)
mean(sim_data$p.obs)
xtabs(~w, data=sim_data)

# condition on being observed
sim_data_obs <- sim_data[ sim_data$w >0,]
dim(sim_data_obs)

nind = nrow(sim_data_obs)
nz   = 1000  # augmented population

data.list <- list(nind=nind, 
                  nz=nz, T=2, 
                  wt=c(sim_data_obs$x, rep(NA, nz)),
                  z =c(rep(1,nind),    rep(NA, nz)),
                  #Y= rbind(as.matrix(sim_data_obs[,c("y1","y2")]), matrix(0, ncol=2, nrow=nz) ,
                  obs= as.numeric(c(sim_data_obs[,c("w"),drop=TRUE]>0, rep(0, nz)))
)

# check the list
data.list




# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

init.list <- list(
      list(a0=1, a1=-1, mu.wt=0, tau.wt=1, psi=.3, wt=c(rep(NA,nind),rnorm(nz,mu_x, sigma_x))),
      list(a0=1, a1=-1, mu.wt=0, tau.wt=1, psi=.3, wt=c(rep(NA,nind),rnorm(nz,mu_x, sigma_x))),
      list(a0=1, a1=-1, mu.wt=0, tau.wt=1, psi=.3, wt=c(rep(NA,nind),rnorm(nz,mu_x, sigma_x)))
)
 

# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("N","a0","a1", "mu.wt", "tau.wt", "sigma", "psi") # parameters to monitor
 

   


# Finally, the actual call to JAGS
set.seed(234234)  # intitalize seed for MCMC 

results <- jags( 
      data      =data.list,   # list of data variables
      inits     =init.list,   # list/function for initial values
      parameters=monitor.list,# list of parameters to monitor
      model.file="model.txt",  # file with bugs model
      n.chains=3,
      n.iter  =50000,          # total iterations INCLUDING burn in
      n.burnin=2000,          # number of burning iterations
      n.thin=20,               # how much to thin
      DIC=TRUE,               # is DIC to be computed?
      working.dir=getwd()    # store results in current working directory
      )


# now results is a BIG list of stuff
names(results)
names(results$BUGSoutput)



#######################################
# extract some of the usual stuff and use R code directly
# use the standard print method
results

# get the summary table
results$BUGSoutput$summary
results$BUGSoutput$summary[,c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]

