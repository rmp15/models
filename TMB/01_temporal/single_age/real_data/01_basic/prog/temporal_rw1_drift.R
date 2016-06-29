library(TMB)

# real data
dat <- readRDS('../data/datus_state_rates_1982_2010')
dat <- subset(dat, fips == 1 & month == 1)
# filter for a single age
age.selected <- 75
sex.selected <- 1
dat <- subset(dat, age==age.selected & sex==sex.selected)

# series length (T) and number of different series (N)
T <- length(unique(dat$year)) * length(unique(dat$month))
N <- length(unique(dat$fips))

# create matrix with data
deaths <- t(matrix(dat$deaths,T,N))
E <- t(matrix(dat$pop.adj,T,N))

# compile cpp file
compile('temporal_rw1_drift.cpp')
dyn.load(dynlib('temporal_rw1_drift'))

# prepare list of parameters for TMB
data <- list(deaths = deaths, E=E)
parameters <- list(beta_0=1.,log_mu=matrix(1.,N,T), pi=matrix(0.,N,T) ,log_prec_rw=1.,log_prec_epsilon=1.)

# run TMB model on simulated data
obj <- MakeADFun(data, parameters, random = c('log_mu','pi'), DLL='temporal_rw1_drift')
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
log_mu <- extract.variable(sd, 'log_mu', 'random')
rw1.pred <- extract.variable(sd,'pi', 'random')
alpha_0 <- summary(sd)[grep("alpha_0", rownames(summary(sd))), 1]
rw1.pred$estimate <- rw1.pred$estimate - alpha_0 # accounts for intercept included in rw1

# obtain joint precision matrix
joint.prec <- sd$jointPrecision

# compare results with INLA
library(INLA)
dat$t <- rep(seq(1,T),N)
dat$t3 <- dat$t2 <- dat$t
dat$e <- seq(1:nrow(dat))
fml <- deaths ~ 1 + 
	f(t, model = "linear", mean.linear = 0, prec.linear = 1e-2) +
	f(t2,model='rw1',hyper = list(theta = list(prior = "loggamma", param = c(1 , 1e-3))))+ 
	f(e, model = 'iid',hyper = list(theta = list(prior = "loggamma", param = c(1 , 1e-3))))

system.time(inla.fit <- inla(fml,family='poisson',data=dat, 	
		control.fixed = list(mean.intercept = 0, prec.intercept = 1e-2),
		E=pop.adj, control.predictor = list(link = 1)))

plot.inla <- inla.fit$summary.fitted.values
plot.inla$id <- seq(1:nrow(plot.inla))

tmbres <- summary(sd)[grep("^mu", rownames(summary(sd))), 1]

inlares <-  inla.fit$summary.fitted.values$`0.5quant`
plot(tmbres, inlares)
abline(a=1,b=1)

# plot fitted values from INLA and TMB compared with original data
library(ggplot2)
p1 <-   ggplot() +
        geom_line(data=dat,colour='blue',aes(x=t,y=rate.adj)) +
        geom_line(data=plot.inla,colour='red',aes(x=id, y=mean)) +
        geom_ribbon(data=plot.inla,alpha=0.3,fill='red',aes(x=id,ymax=(`0.975quant`),ymin=(`0.025quant`))) +
        geom_line(data=log_mu,colour='green',aes(x=id,y=exp(estimate))) +
        geom_ribbon(data=log_mu,alpha=0.2,fill='green',
			aes(x=id,ymax=exp(estimate+1.96*std.error),ymin=exp(estimate-1.96*std.error))) +
        ggtitle('Fitted values, blue=data, red=INLA, green=TMB') +
        theme_bw()

# plot random walks from INLA and TMB
plot.inla.rw1 <- inla.fit$summary.random$t2
plot.inla.rw1$id <- seq(1:nrow(plot.inla.rw1))

p2 <-   ggplot() +
        geom_line(data=plot.inla.rw1,colour='red',aes(x=id,y=mean)) +
        geom_ribbon(data=plot.inla.rw1,alpha=0.2,fill='red',aes(x=id,ymax=(`0.975quant`),ymin=(`0.025quant`))) +
        geom_line(data=rw1.pred,colour='green',aes(x=id,y=estimate)) +
        geom_ribbon(data=rw1.pred,fill='green',alpha=0.2,aes(x=id,ymax=estimate+1.96*std.error,ymin=estimate-1.96*std.error)) +
        ggtitle('Random Walk: fitted values, red=INLA, green=TMB') +
        theme_bw()

pdf('results.pdf',paper='a4r',height=0,width=0)
p1
p2
dev.off()

file.remove("temporal_rw1_drift.o", "temporal_rw1_drift.so")
