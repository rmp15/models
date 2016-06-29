library(TMB)
library(dplyr)

# real data
dat <- readRDS('../data/datus_state_rates_1982_2010')
# filter for a single age at a single time
age.selected <- 65
sex.selected <- 1
year.selected <- 1982
month.selected <- c(1:12)
dat <- subset(dat, age==age.selected & sex==sex.selected & 
		year==year.selected & month==month.selected)

# extract unique table of year and months to generate year.month
dat.year.month <- unique(dat[,c('year', 'month')])
dat.year.month$month <- as.integer(dat.year.month$month)
dat.year.month$year.month <- seq(nrow(dat.year.month))

# merge year.month table with population table to create year.month id
dat <- merge(dat,dat.year.month, by=c('year','month'))

# series length (T) and number of different series (N)
T <- length(unique(dat$year)) * length(unique(dat$month))
N <- length(unique(dat$fips))

# load maptools and USA map
library(maptools)
library(RColorBrewer)
getinfo.shape("../data/shapefiles/states")
USA.gen <- readShapePoly("../data/shapefiles/states")
#plot(USA.gen)

# extract data from shapefile
shapefile.data <- attr(USA.gen, 'data')
names(shapefile.data)[3] <- 'fips'
shapefile.data$fips <- as.integer(as.character(shapefile.data$fips))

# re-insert back into shapefile
attr(USA.gen,'data') <- shapefile.data

# create lookup for fips and DRAWSEQ
drawseq.lookup <- as.data.frame(cbind(DRAWSEQ=shapefile.data$DRAWSEQ,fips=shapefile.data$fips))

# load rgdal and spdep
#library(rgdal)
library(spdep)

# create adjacency matrix
USA.nb <- poly2nb(USA.gen, queen=1)
#plot(USA.nb,coordinates(USA.gen),add=1)

# make matrix compatible with INLA
library(INLA)

nb2INLA("USA.graph",USA.nb)
USA.adj <- "USA.graph"

# Add connections Hawaii -> California, Alaska -> Washington
USA.adj <- "USA.graph.edit"

# create graph from shapefile, plot if desired
H <- inla.read.graph(filename=USA.adj)
#image(inla.graph2matrix(H),xlab="",ylab="")

# create sparse matrix and make diagonals zero
USA.adj.matrix <- inla.graph2matrix(H)
diag(USA.adj.matrix) <- 0

# create diagonal matrix with number of neighbours on diagonal
USA.diag.matrix <- diag(rowSums(USA.adj.matrix))

# create precision matrix
prec.mat <- USA.diag.matrix - USA.adj.matrix
prec.mat <- as(prec.mat, "dgTMatrix")

dat <- merge(dat,drawseq.lookup, by='fips')
dat <- dat[order(dat$DRAWSEQ,dat$sex,dat$age,dat$year.month),]

# create matrix with data
deaths <- t(matrix(dat$deaths,T,N))
E <- t(matrix(dat$pop.adj,T,N))
P <- as(prec.mat, 'dgTMatrix')

# compile cpp file
TMB::compile("spatiotemporal.cpp")
dyn.load(dynlib("spatiotemporal"))

# prepare list of parameters for TMB
data <- list(deaths = deaths, E=E, P=P)
parameters <- list(alpha = 0., V = rep(0., N), W = rep(0., N), log_sigma2_V = 0., log_sigma2_U = 0.)

obj <- MakeADFun(
	data = data,
	parameters = parameters, 
	random = c("V", "W"),
	DLL = "spatiotemporal",
	hessian = TRUE
)
obj$hessian <- TRUE
opt <- do.call('optim', obj)
rep <- sdreport(obj)
tmb.res <- summary(rep)

# extract predicted values from TMB
deaths_pred <- tmb.res[rownames(tmb.res) == "deaths_pred", "Estimate"]
deaths_pred.sd <- tmb.res[rownames(tmb.res) == "deaths_pred", "Std. Error"]
mu.tmb <- tmb.res[rownames(tmb.res) == "mu", "Estimate"]
mu.tmb.sd <- tmb.res[rownames(tmb.res) == "mu", "Std. Error"]

# Fit in INLA
library(INLA)
dat.inla <- data.frame(dat)
#dat.inla <- merge(dat.inla,drawseq.lookup, by='fips')

# add columns for INLA
dat.inla$ID2 <- dat.inla$ID <- dat.inla$DRAWSEQ
dat.inla$year.month2 <- dat.inla$year.month
dat.inla$e <- seq(1:nrow(dat.inla))

fml <- deaths ~ 1 +
        	year.month +                                               
		f(ID, model = "iid", hyper = list(theta = 
		list(prior = "loggamma", param = c(0.5, 5e-4)))) + 
		f(ID2, model = "besag", graph = USA.adj,
		hyper = list(theta = list(prior = "loggamma", 
		param = c(0.5 , 5e-4)))) +
        	f(year.month2, model="rw1") + 	
		f(e, model = 'iid',hyper = list(theta = list(prior = "loggamma", 
		param = c(1 , 1e-3))))

inla.res <- inla(
	formula=fml,
	data = dat.inla,
	family = "poisson", 
	E = pop.adj,
	control.predictor = list(link = 1),
	control.compute = list(dic = TRUE, config = TRUE),
	control.fixed = list(mean.intercept = 0, prec.intercept = 1e-2), 
)

# combine results into a dataset
dat <- data.frame(dat, rate.inla.mean=inla.res$summary.fitted.values$mean, 
		rate.inla.sd=inla.res$summary.fitted.values$sd,
		rate.tmb.mean=mu.tmb,
		rate.tmb.sd=mu.tmb.sd)

library(ggplot2)
ggplot(data=dat,aes(x=DRAWSEQ)) +
geom_line(aes(y=rate.inla.mean),color='red') +
geom_ribbon(fill='red',alpha=0.2,aes(ymax=rate.inla.mean+1.96*rate.inla.sd,ymin=rate.inla.mean-1.96*rate.inla.sd)) +
geom_line(aes(y=rate.tmb.mean),color='green') +
geom_ribbon(fill='green',alpha=0.2,aes(ymax=rate.tmb.mean+1.96*rate.tmb.sd,ymin=rate.tmb.mean-1.96*rate.tmb.sd)) +
xlab('state') +
ylab('death rate') +
theme_bw()

# TMB vs INLA estimates
# Intercept
tmb.res["alpha", "Estimate"] 
inla.res$summary.fixed["(Intercept)", "mean"]
# Variance of V
exp(tmb.res["log_sigma2_V", "Estimate"])
1 / inla.res$summary.hyperpar["Precision for ID", "mean"]
# Variance of U
exp(tmb.res["log_sigma2_U", "Estimate"])
1 / inla.res$summary.hyperpar["Precision for ID2", "mean"]

# fitted values
plot(
	dat$rate, 
	mu.tmb,
	cex = 0.5, 
	lwd = 0.3,
	col = "red", 
	xlab = "unsmoothed values", 
	ylab = "TMB (red) / INLA (blue)"
)
points(
	dat$rate,
	inla.res$summary.fitted.values$mean,
	cex = 0.5, 
	lwd = 0.3,
	col = "blue"
)
abline(0, 1, col = "grey", lty = 2)

# TMB vs INLA
plot(
	inla.res$summary.fitted.values$mean, 
	mu.tmb,
	cex = 0.5, 
	lwd = 0.3,
	col = "red", 
	xlab = "unsmoothed values", 
	ylab = "TMB (red) / INLA (blue)"
)
abline(0, 1, col = "grey", lty = 2)

# U - Spatial term
plot(
	tmb.res[rownames(tmb.res) == "U", "Estimate"],
	inla.res$summary.random$ID2$mean,
	cex = 0.5, 
	lwd = 0.3,
	col = "blue"
)
abline(0, 1, col = "grey", lty = 2)

# V - Non-spatial term
plot(
	tmb.res[rownames(tmb.res) == "V", "Estimate"],
	inla.res$summary.random$ID$mean,
	cex = 0.5, 
	lwd = 0.3,
	col = "blue"
)
abline(0, 1, col = "grey", lty = 2)

file.remove("spatiotemporal.o", "spatiotemporal.so")

