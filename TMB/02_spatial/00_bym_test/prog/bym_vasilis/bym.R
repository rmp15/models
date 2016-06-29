# Generate some data
adj <- eval(parse(file = "adjacency.txt"))
adj.list <- Map(function(i1, i2) adj$adj[i1:i2], 
	i1 = cumsum(adj$num)-adj$num + 1,
	i2 = cumsum(adj$num))
AI <- do.call(cbind,
	lapply(adj.list, function(x) as.numeric(seq_along(adj$num) %in% x)))
dat <- eval(parse(file = "lung_data.txt"))

DI <-  diag(rowSums(AI))
prec.mat <- DI - AI

# Fit model in TMB
library(TMB)
compile("bym.cpp")
dyn.load(dynlib("bym"))
obj <- MakeADFun(
	data = list(E = dat$E, O = dat$O, P = as(prec.mat, "dgTMatrix")),
	parameters = list(alpha = 0., V = rep(0., dat$N), W = rep(0., dat$N), log_sigma2_V = 0., log_sigma2_U = 0.), 
	random = c("V", "W"),
	DLL = "bym",
	hessian = TRUE
)
obj$hessian <- TRUE
opt <- do.call('optim', obj)
rep <- sdreport(obj)
tmb.res <- summary(rep)

# Fit in INLA
library(INLA)
inla.dat <- data.frame(dat)
inla.dat$i <- inla.dat$j <- seq(nrow(inla.dat))
inla.res <- inla(
	O ~ f(j, model = "iid", hyper = list(theta = list(prior = "loggamma", param = c(0.5, 5e-4)))) + 
		f(i, model = "besag", graph = inla.graph2matrix(AI),
			hyper = list(theta = list(prior = "loggamma", param = c(0.5 , 5e-4)))),
	data = inla.dat,
	family = "poisson", 
	E = E,
	control.predictor = list(link = 1),
	control.compute = list(dic = TRUE, config = TRUE),
	control.fixed = list(mean.intercept = 0, prec.intercept = 1e-2), 
)

# OpenBUGS results
ob.res <- readRDS("OB_BYM.rds")

# TMB vs OpenBUGS vs INLA estimates
# Intercept
tmb.res["alpha", "Estimate"] 
inla.res$summary.fixed["(Intercept)", "mean"]
ob.res["alpha", "mean"]
# Variance of V
exp(tmb.res["log_sigma2_V", "Estimate"])
ob.res["sigma2.v", "mean"]
1 / inla.res$summary.hyperpar["Precision for j", "mean"]
# Variance of U
exp(tmb.res["log_sigma2_U", "Estimate"])
ob.res["sigma2.u", "mean"]
1 / inla.res$summary.hyperpar["Precision for i", "mean"]

# RR - fitted values
plot(
	tmb.res[rownames(tmb.res) == "RR", "Estimate"],
	ob.res[grep("^RR", rownames(ob.res)), "mean"], 
	cex = 0.5, 
	lwd = 0.3,
	col = "red", 
	xlab = "TMB", 
	ylab = "OpenBUGS (red) / INLA (blue)"
)
points(
	tmb.res[rownames(tmb.res) == "RR", "Estimate"],
	inla.res$summary.fitted.values$mean,
	cex = 0.5, 
	lwd = 0.3,
	col = "blue"
)
abline(0, 1, col = "grey", lty = 2)

# U - Spatial term
plot(
	tmb.res[rownames(tmb.res) == "U", "Estimate"],
	ob.res[grep("^U\\[", rownames(ob.res)), "mean"], 
	cex = 0.5, 
	lwd = 0.3,
	col = "red", 
	xlab = "TMB", 
	ylab = "OpenBUGS (red) / INLA (blue)"
)
points(
	tmb.res[rownames(tmb.res) == "U", "Estimate"],
	inla.res$summary.random$i$mean,
	cex = 0.5, 
	lwd = 0.3,
	col = "blue"
)
abline(0, 1, col = "grey", lty = 2)

# V - Non-spatial term
plot(
	tmb.res[rownames(tmb.res) == "V", "Estimate"],
	ob.res[grep("^V\\[", rownames(ob.res)), "mean"], 
	cex = 0.5, 
	lwd = 0.3,
	col = "red", 
	xlab = "TMB", 
	ylab = "OpenBUGS (red) / INLA (blue)"
)
points(
	tmb.res[rownames(tmb.res) == "V", "Estimate"],
	inla.res$summary.random$j$mean,
	cex = 0.5, 
	lwd = 0.3,
	col = "blue"
)
abline(0, 1, col = "grey", lty = 2)

file.remove("bym.o", "bym.so")
