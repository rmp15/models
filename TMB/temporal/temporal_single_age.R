rm(list=ls())

set.seed(123)

library(INLA)
library(TMB)

# create a dataset from specified model

# series length (t) and number of different series (n)
T <- 10000
N <- 1

# parameters
alpha_0 <- 10
beta_0 <- -0.01
prec_rw1 <- 150000

# generate random walk
generate.rw1 <- function(prec) {
dummy <-c(0,cumsum(rnorm(n=T-1, mean=0,sd=1/sqrt((prec)))))
return(dummy)
}

rw1 <- replicate(N,generate.rw1(prec_rw1)) 

# generate covariate values
t <- rep(seq(1:T),N)

# compute mu values
real_mu <- exp(alpha_0 + beta_0 * t  + rw1)

y <- rpois(n=T*N, lambda=real_mu)

dat <- data.frame(t=t,mu=y)

# fit using glm
# glm.fit <- glm(dat$mu ~ 1+t, family='poisson'(link = "log"))

# fit using INLA
t2 <- t
fml <- mu ~ 1 + t + f(t2,model='rw1')
inla.fit <- inla(fml,family='poisson',data=dat)

# create matrix with data
log_mu <- matrix(dat$log_mu,t,n)

