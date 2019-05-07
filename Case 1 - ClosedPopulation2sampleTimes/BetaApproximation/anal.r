# Compare the E[1/p] vs 1/E[p] for a beta distribution

library(ggplot2)


gen.ratio <- function(alpha,beta,n.sim=10000){

   dist<- data.frame(p=rbeta(n.sim, alpha, beta))

   plot1 <- ggplot(dist, aes(x=p))+
     geom_histogram(alpha=0.2)
   plot(plot1)

  # theoretical mean
  t.mean = alpha/(alpha+beta)
 
  # empirical mean
  e.mean <- mean(dist$p)

  # theoretical inv of mean
  t.inv.mean <- 1/t.mean

  print(summary(dist$p))
  print(summary(1/dist$p))

  plot2 <-ggplot(dist, aes(x=1/p))+
    geom_histogram(alpha=0.2)
  plot(plot2)

  # empirical e[1/p]
  e.mean.inv.p <- mean(1/dist$p)

  # theoretical t[1/p]
  t.mean.inv.p <- (alpha+beta-1)/(alpha-1)

  # compute relative increase in variance
  rel.increase.var <-( t.mean * t.mean.inv.p - 1 )/(1-t.mean)

  data.frame(alpha, beta, t.mean, e.mean, t.inv.mean, t.mean.inv.p, e.mean.inv.p, rel.increase.var )
}


###################################################################################
###################################################################################
###################################################################################

gen.ratio(2, 8)

gen.ratio(1.1, 8.9)

gen.ratio(2,2)
###################################################################################
###################################################################################
###################################################################################

gen.ratio(10,40)

###################################################################################
###################################################################################
###################################################################################

gen.ratio(20,80)

###################################################################################
###################################################################################
###################################################################################

gen.ratio(200,800)


###################################################################################
###################################################################################
###################################################################################
# what happens in Case 1

dist <- data.frame(x=rnorm(100000,0,1))
dist$logitp <- -1 + 1*dist$x
dist$p1     <- 1/(1+exp(-dist$logitp))
dist$p2     <- dist$p1
dist$p      <- 1-(1-dist$p1)*(1-dist$p2)
ggplot(data=dist, aes(x=p1))+
  ggtitle ("Distribution of capture prob on each occasion")+
  geom_histogram(alpha=0.2)


plot1 <- ggplot(data=dist, aes(x=p, y=..density..))+
   ggtitle("Distribution of overall prob of capture over both occasions")+
   geom_histogram(alpha=0.2)
plot1

mean(dist$p)
mean(1/dist$p)

# fit a beta distribution to p
library(fitdistrplus)
parms <- fitdist(dist$p, "beta")
parms$estimate
fit.beta <- data.frame(p=seq(0,1,.01), dp=dbeta(seq(0,1,.01), parms$estimate["shape1"], parms$estimate["shape2"]))
plot1 + geom_line(data=fit.beta, aes(y=dp))


parms$estimate["shape1"]/sum(parms$estimate)
(sum(parms$estimate)-1)/(parms$estimate[1]-1)

(mean(dist$p)*mean(1/dist$p)-1)/(1-mean(dist$p))

gen.ratio(parms$estimate[1], parms$estimate[2])

gen.ratio(parms$estimate[1], parms$estimate[2])
gen.ratio(1.85, 1.72)
gen.ratio(1.76, 1.84)

# Huggins 1048, 95% ci 857 to 1342 -> se= 123
# JAGS    1073  95% ci SD 134
# MLE     1052  SE  125
# Petersen 1115 SE  63

(134^2-125^2)/125^2

# petersen estimate se=63
# jags estimate sd=136

(136/63)^2
