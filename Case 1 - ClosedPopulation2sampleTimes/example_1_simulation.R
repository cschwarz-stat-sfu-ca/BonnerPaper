## Load Packages
library(tidyverse)
library(rjags)
library(RMark)

## Set parameters

N <- 10000 # Population size
mu_x <- 0                               # Covariate mean
sigma_x  <- 1                           # Covaraite sd

beta <- c(-1,1)                          # Parameters of capture probability
#beta <- c(0,0)
## Set seed
set.seed(8888)

## Simulate data
sim_data <- tibble(x = rnorm(N,mu_x,sigma_x),
                   p = (1 + exp(-(beta[1] + beta[2]*x)))^-1,
                   y1 = rbinom(N,1,p),
                   y2 = rbinom(N,1,p),
                   w = y1 + y2,
                   p.obs= 1-(1-p)^2)
mean(sim_data$p.obs)
xtabs(~w, data=sim_data)

# get the simple petersen estimator
library(BTSPAS)
SimplePetersen(n1=sum(sim_data$y1), u2=sum((1-sim_data$y1)*sim_data$y2), m2=sum(sim_data$y1*sim_data$y2))

1/mean(sim_data$p.obs)
mean(1/sim_data$p.obs)

## Extract data for detected individuals
obs_data <- filter(sim_data,w > 0)


# Do a full maximum likelihood
# parameters <- beta0, beta1, mu, sigma, N
# we will use logsigma and logN as parameter

library(bbmle)

loglik <- function(beta0, beta1, mu, logsigma, logN, history){
   # actual log likelihood for case 1 of Simon's work
   # par = parameters to optmize beta0, beta1, mu, logsigma, logN
   sigma <- exp(logsigma)
   N     <- exp(logN)
   n     <- nrow(history)
   # compute p00 = probability that an animal not captured
   pnotcapture <- function(x, beta0, beta1, mu, sigma){
       # logit capture probability
       logitp <- beta0 + beta1*x
       p <- 1/(1+exp(-logitp))
       res <- (1-p)*(1-p)
       res <- res*dnorm(x, mean=mu, sd=sigma)
       #browser()
       res
   }
   p00 <- integrate(pnotcapture,  lower=-100, upper=100, beta0=beta0, beta1=beta1, mu=mu, sigma=sigma)

   #browser()
   # Contribution from captured animals
   logitp <- beta0 + beta1*history$x
   p <- 1/(1+exp(-logitp))
   logphist <-history$y1*log(p) + (1-history$y1)*log(1-p) + history$y2*log(p) + (1-history$y2)*log(1-p)+
              dnorm(history$x, mean=mu, sd=sigma, log=TRUE)
   
   # contribution from uncaptured animals
   logp00  <- (N-n)*log(p00$value)
   
   # contribution from factorials
   logfact <- lgamma(N+1)-lgamma(N-n+1)
   
   # add everything up
   loglik <- sum(logphist)+logp00 + logfact
   
   # return negative  log likelihood
   #cat (beta0, beta1, mu, sigma, N, n, p00$value, sum(logphist), logp00, logfact, loglik, "\n")

   -loglik
}



loglik(-1, 1, 0, log(sigma_x), log(N), obs_data[,c("y1","y2","x")])


fit <- mle2( loglik, start=list(beta0=beta[1], beta1=beta[2], mu=mu_x, logsigma=log(sigma_x), logN=log(N)),
          #fixed=list(mu=0), #control=list(maxit=5),
          #lower=c(-1,-1,-1,-1,6.5),upper=c(1,1,1,1,10),
          data=list(history=obs_data[,c("y1","y2","x")]))
summary(fit)
coef(fit)
stdEr(fit)
exp(coef(fit)["logN"])
exp(coef(fit)["logN"])*stdEr(fit)["logN"]

vcov(fit)
cov2cor(vcov(fit))




ggplot(data=sim_data, aes(x=p.obs))+
   ggtitle("Distribution of p(observing animal)")+
   geom_histogram(alpha=0.2)
library(fitdistrplus)
fitdist(sim_data$p.obs, "beta")



## Run model in jags
jags_data <- list(n=nrow(obs_data),
                  w=pull(obs_data,"w"), # Number of captures
                  x=pull(obs_data,"x"), # Covariate
                  dummy=rep(0,nrow(obs_data)))

jags_inits <- list(beta=c(0,0))

jags_model_1 <- jags.model("example_1_model.jag",
                           jags_data,
                           jags_inits,
                           n.adapt=1000)

jags_samples_1  <- coda.samples(jags_model_1,
                                c("beta","N"),
                                5000)

summary(jags_samples_1)

## Run Huggins model in MARK for comparison
mark_data <- data.frame(ch=paste0(obs_data$y1,obs_data$y2),
                        x=obs_data$x,
                        freq=rep(1,nrow(obs_data)),
                        stringsAsFactors = FALSE)

mark_proc <- process.data(mark_data,model="Huggins")

mark_ddl <- make.design.data(mark_proc)

p_x <- list(formula=~x,share=TRUE)

mark_fit <- mark(mark_proc,mark_ddl,model.parameters = list(p=p_x))
mark_fit$results


# Huggins 1048, 95% ci 857 to 1342 -> se= 123
# JAGS    1073  95% ci SD 134
# MLE     1052  SE  125
# Petersen 1115 SE  63



