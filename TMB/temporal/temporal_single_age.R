rm(list=ls())

library(INLA)
library(TMB)

# create a dataset from specified model

# series length (t) and number of different series (n)
t <- 100
n <- 1

# coefficients
alpha_0 <- 1
beta_0 <- 0.01

# generate random walk
generate.rw1 <- function(log_sigma_rw) {
dummy <-c(0,cumsum(rnorm(n=t-1, mean=0,sd=exp(log_sigma_rw))))
return(dummy)
}

rw1 <- rep(generate.rw1(0),n) 

# generate covariate values
x <- rep(seq(1:t),n)

# compute mu values
real_mu <- exp(alpha_0 + beta_0 * x + rw1)

y <- rpois(n=t, lambda=real_mu)

dat <- data.frame(x=x,mu=y,log_mu=log(y))

# fit using built-in glm
fit <- glm(mu~x, family= poisson, data=dat)

# create matrix with data


