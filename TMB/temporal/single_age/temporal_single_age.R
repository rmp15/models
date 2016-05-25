library(TMB)

# series length (t) and number of different series (n)
T <- 5000
N <- 1

# parameters
alpha_0 <- 10
beta_0 <- -0.001
prec_rw1 <- 100
prec_e <- 300

# generate random walk
generate.rw1 <- function(prec) {
#set.seed(123)
dummy <-c(0,cumsum(rnorm(n=T-1, mean=0,sd=1/sqrt((prec)))))
return(dummy)
}

# generate random walk and overdispersion
rw1 <- replicate(N,generate.rw1(prec_rw1)) 
od  <- rnorm(T, mean = 0, sd = sqrt(1 / prec_e))

# generate covariate values
t <- rep(seq(1:T),N)

# compute mu values
real_lambda <- exp(alpha_0 + beta_0 * t  + rw1 + od)

counts <- rpois(n=T*N, lambda=real_lambda)

dat <- data.frame(t=t,counts=counts)

# create matrix with data
log_counts <- matrix(log(dat$counts),N,T)

# compile cpp file
compile('model.cpp')
dyn.load(dynlib('model'))

# prepare list of parameters for TMB
data <- list(log_counts = log_counts)
parameters <- list(beta_0=1.,log_tau_rw=4.,log_tau_epsilon=1.,log_counts_pred=matrix(1.,N,T), pi=matrix(0.,N,T) )

# run TMB model on simulated data
#obj <- MakeADFun(data, parameters, random = "pi", DLL='model')
obj <- MakeADFun(data, parameters, random = "pi", DLL='model')
obj$hessian <- FALSE
opt <- do.call("optim", obj)
sd <- sdreport(obj)

# extract fixed effects
fixed <- summary(sd, 'fixed')

# compare results with INLA
library(INLA)
t2 <- t
fml <- counts ~ 1 + t + f(t2,model='rw1') + f(t3, model = "iid")
inla.fit <- inla(fml,family='poisson',data=dat)
