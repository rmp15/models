library(TMB)

# series length (t) and number of different series (n)
T <- 100
N <- 1

# parameters
alpha_0 <- 10
beta_0 <- -0.01

# generate covariate values
t <- rep(seq(1:T),N)

# compute mu values
real_lambda <- exp(alpha_0 + beta_0 * t)

counts <- rpois(n=T*N, lambda=real_lambda)

dat <- data.frame(t=t,counts=counts)

# create matrix with data
log_counts <- matrix(log(dat$counts),N,T)

# compile cpp file
compile('temporal_drift.cpp')
dyn.load(dynlib('temporal_drift'))

# prepare list of parameters for TMB
data <- list(log_counts = log_counts)
parameters <- list(alpha,=1.,beta_0=1.,log_counts_pred=matrix(1.,N,T))

# run TMB model on simulated data
obj <- MakeADFun(data, parameters, DLL='temporal_drift')
obj$hessian <- FALSE
opt <- do.call("optim", obj)
sd <- sdreport(obj)

# extract fixed effects
fixed <- summary(sd, 'fixed')

# compare results with INLA
library(INLA)
fml <- counts ~ 1 + t
inla.fit <- inla(fml,family='poisson',data=dat, control.predictor = list(link = 1))

# plot results compared with original data
library(ggplot2)
plot.inla <- inla.fit$summary.fitted.values
plot.inla$id <- seq(1:nrow(plot.inla))
p <-  ggplot() +
      geom_line(data=dat,colour='blue',aes(x=t,y=counts)) + 
      geom_line(data=plot.inla,colour='red',aes(x=id, y=mean)) + 
      geom_ribbon(data=plot.inla,alpha=0.5,aes(x=id,ymax=(`0.975quant`),ymin=(`0.025quant`)))
