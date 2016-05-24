rm(list=ls())

library(INLA)
library(TMB)

# create a dataset from specified model

# series length (t) and number of different series (n)
T <- 500
N <- 1

# parameters
alpha_0 <- 10
beta_0 <- -0.01
prec_rw1 <- 150000

# generate random walk
generate.rw1 <- function(prec) {
set.seed(123)
dummy <-c(0,cumsum(rnorm(n=T-1, mean=0,sd=1/sqrt((prec)))))
return(dummy)
}

rw1 <- replicate(N,generate.rw1(prec_rw1)) 

# generate covariate values
t <- rep(seq(1:T),N)

# compute mu values
real_lambda <- exp(alpha_0 + beta_0 * t  + rw1)

counts <- rpois(n=T*N, lambda=real_lambda)

dat <- data.frame(t=t,counts=counts)

# fit using glm
# glm.fit <- glm(dat$y ~ 1+t, family='poisson'(link = "log"))

# fit using INLA
t2 <- t
fml <- counts ~ 1 + t + f(t2,model='rw1')
inla.fit <- inla(fml,family='poisson',data=dat)

# create matrix with data
log_counts <- matrix(log(dat$counts),T,N)

# compile cpp file
compile('temporal_single_age.cpp')
dyn.load(dynlib('temporal_single_age'))

data <- list(log_counts = log_counts)
parameters <- list(alpha_0=0,beta_0=0,log_tau_rw=-1,log_tau_epsilon=-1)

obj <- MakeADFun(data, parameters,random='log_counts_pred', DLL='temporal_single_age')
