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
# glm.fit <- glm(dat$counts ~ 1+t, family='poisson'(link = "log"))

# fit using INLA
t2 <- t
fml <- counts ~ 1 + t + f(t2,model='rw1')
inla.fit <- inla(fml,family='poisson',data=dat)

# create matrix with data
log_counts <- matrix(log(dat$counts),N,T)

# compile cpp file
compile('temporal_single_age.cpp')
dyn.load(dynlib('temporal_single_age'))

# prepare list of parameters for TMB
data <- list(log_counts = log_counts)
parameters <- list(alpha_0=10,beta_0=-0.02,log_tau_rw=-1,log_tau_epsilon=-1,log_counts_pred=matrix(0,N,T))

# run TMB model on simulated data
obj <- MakeADFun(data, parameters, DLL='temporal_single_age')
obj$hessian <- TRUE
opt <- do.call("optim", obj)
sd <- sdreport(obj)

# extract fixed effects
fixed <- summary(sd, 'fixed')
