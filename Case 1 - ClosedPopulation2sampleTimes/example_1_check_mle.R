# Check my coding of the mle via simulation ...

## Load Packages
library(tidyverse)
library(rjags)
library(RMark)
library(bbmle)

## Set parameters

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


sim_mle <- function(N, mu_x, sigma_x, beta){
  
  # simulate some data
  sim_data <- tibble(x = rnorm(N,mu_x,sigma_x),
                   p = (1 + exp(-(beta[1] + beta[2]*x)))^-1,
                   y1 = rbinom(N,1,p),
                   y2 = rbinom(N,1,p),
                   w = y1 + y2,
                   p.obs= 1-(1-p)^2)
  
   ## Extract data for detected individuals
   obs_data <- filter(sim_data,w > 0)
 
   ## do the fit
   fit <- tryCatch({mle2( loglik, start=list(beta0=beta[1], beta1=beta[2], mu=mu_x, logsigma=log(sigma_x), logN=log(N)),
          #fixed=list(mu=0), #control=list(maxit=5),
          #lower=c(-1,-1,-1,-1,6.5),upper=c(1,1,1,1,10),
          data=list(history=obs_data[,c("y1","y2","x")]))},
          
          error=function(cond) {
            message("MLE failed to converge\n")
            message(cond)
            # Choose a return value in case of error
            return(NULL)
        }
    )    
   ## return the results
   list(N=N, mu_x=mu_x, sigma_x=sigma_x, beta=beta, fit=fit)
  }

set.seed(8888)

all.fits <- plyr::llply(1:100, function(sim){
   fit <- sim_mle(N=1000, mu_x=0, sigma_x=1, beta=c(-1,1))
   fit$sim = sim
   fit
})

# extract the estimates from the fits
est.summary <- plyr::ldply(all.fits, function(x){
   cat("extracting from ", x$sim, "\n")
   #if(x$sim ==8) browser()
   if(!is.null(x$fit)){
        res <- data.frame(Nhat=exp(coef(x$fit)["logN"]),
              Nhat.se=exp(coef(x$fit)["logN"])*stdEr(x$fit)["logN"])
   }
   if(is.null(x$fit)){
        res<-data.frame(Nhat=NA,
              Nhat.se=NA)

   }
   res
})
est.summary

est.summary <- est.summary[complete.cases(est.summary),]

apply(est.summary,2,mean)
apply(est.summary,2,sd)
hist(est.summary$Nhat)
