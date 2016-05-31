library(TMB)

# series length (T) and number of different series (N)
T <- 10000
N <- 1

# parameters
alpha_0 <- 15
beta_0 <- -0.001
prec_rw1 <- 100
prec_e <- 300

# generate random walk
generate.rw1 <- function(prec) {
set.seed(123)
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
compile('temporal_rw1_drift.cpp')
dyn.load(dynlib('temporal_rw1_drift'))

# prepare list of parameters for TMB
data <- list(log_counts = log_counts)
parameters <- list(beta_0=1.,log_counts_pred=matrix(1.,N,T), pi=matrix(0.,N,T) ,log_prec_rw=1.,log_prec_epsilon=1.)

# run TMB model on simulated data
obj <- MakeADFun(data, parameters, random = 'pi', DLL='temporal_rw1_drift')
obj$hessian <- FALSE
system.time(opt <- do.call('optim', obj))
system.time(sd <- sdreport(obj,getJointPrecision = TRUE))

# function to extract desired variable
extract.variable <- function(sd, var, type) {
    point <- summary(sd, type)
    point <- point[rownames(point) == var, ][, 'Estimate']
    point <- data.frame(estimate=point)

    error <- summary(sd, type)
    error <- error[rownames(error) == var, ][, 'Std. Error']
    error <- data.frame(std.error=error)
    complete <- cbind(point,error)
    complete$id <- seq(1:nrow(complete))
    return(complete)
}

# extract predicted values from TMB
log_counts.pred <- extract.variable(sd, 'log_counts_pred', 'fixed')
rw1.pred <- extract.variable(sd,'pi', 'random')
rw1.pred$estimate <- rw1.pred$estimate - alpha_0 # accounts for intercept included in rw1

# obtain joint precision matrix
joint.prec <- sd$jointPrecision

# compare results with INLA
library(INLA)
t3 <- t2 <- t
fml <- counts ~ 1 + t + f(t2,model='rw1') + f(t3, model = 'iid')
system.time(inla.fit <- inla(fml,family='poisson',data=dat, control.predictor = list(link = 1)))

plot.inla <- inla.fit$summary.fitted.values
plot.inla$id <- seq(1:nrow(plot.inla))

# plot fitted values from INLA and TMB compared with original data
library(ggplot2)
p1 <-   ggplot() +
        geom_line(data=dat,colour='blue',aes(x=t,y=counts)) +
        geom_line(data=plot.inla,colour='red',aes(x=id, y=mean)) +
        geom_ribbon(data=plot.inla,alpha=0.3,fill='red',aes(x=id,ymax=(`0.975quant`),ymin=(`0.025quant`))) +
        geom_line(data=log_counts.pred,colour='green',aes(x=id,y=exp(estimate))) +
        geom_ribbon(data=log_counts.pred,alpha=0.2,fill='green',aes(x=id,ymax=exp(estimate+1.6449*std.error),ymin=exp(estimate-1.6449*std.error))) +
        ggtitle('Fitted values, red=INLA, green=TMB') +
        theme_bw()

# plot random walks from INLA and TMB compared with original data
rw1.df <- as.data.frame(rw1)
rw1.df$id <- seq(1:nrow(rw1.df))
names(rw1.df) <- c('values','id')

plot.inla.rw1 <- inla.fit$summary.random$t2
plot.inla.rw1$id <- seq(1:nrow(plot.inla.rw1))

p2 <-   ggplot() +
        geom_line(data=rw1.df,aes(x=id,y=values)) +
        geom_line(data=plot.inla.rw1,colour='red',aes(x=id,y=mean)) +
        geom_ribbon(data=plot.inla.rw1,alpha=0.2,fill='red',aes(x=id,ymax=(`0.975quant`),ymin=(`0.025quant`))) +
        geom_line(data=rw1.pred,colour='green',aes(x=id,y=estimate)) +
        geom_ribbon(data=rw1.pred,fill='green',alpha=0.2,aes(x=id,ymax=estimate+1.6449*std.error,ymin=estimate-1.6449*std.error)) +
        ggtitle('Random Walk: fitted values, red=INLA, green=TMB') +
        theme_bw()



