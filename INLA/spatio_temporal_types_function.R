st=ls())

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)

# break down the arguments from Rscript
age.arg <- as.numeric(args[1])
sex.arg <- as.numeric(args[2])
year.start.arg <- as.numeric(args[3])
year.end.arg <- as.numeric(args[4])
type.arg <- as.numeric(args[5])
cluster.arg <- as.numeric(args[6])

require(mailR)

# load US data
dat.inla.load <- readRDS('mortality/datus_state_rates_1982_2010')

# available ages
age.filter <- unique(dat.inla.load$age)

library(dplyr)

# state lookup
state.lookup <- read.csv('fips_lookup/name_fips_lookup.csv')

# gender code
sex.lookup <- c('male','female')

# month lookup
month.lookup <- c('January','February','March','April','May','June','July','August','September','October','November','December')

# load maptools and USA map
library(maptools)
library(RColorBrewer)
getinfo.shape("shapefiles/states")
USA.gen <- readShapePoly("shapefiles/states")
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

# plot matrix if desired
#H <- inla.read.graph(filename=USA.adj)
#image(inla.graph2matrix(H),xlab="",ylab="")

##############

# function to enable age group and sex to be selected
inla.function <- function(age.sel=55,sex.sel=1,year.start=1982,year.end=2010,type=1,cluster=0) {

dat.inla <- dat.inla.load

# filter all data by sex age and month
sex <- sex.sel
age <- age.sel
fit.years <- year.start:year.end
dat.inla <- dat.inla[dat.inla$sex==sex & dat.inla$age==age & dat.inla$year %in% fit.years,]

dat.inla <- merge(dat.inla,drawseq.lookup, by='fips')

# extract unique table of year and months to generate year.month
dat.year.month <- unique(dat.inla[,c('year', 'month')])
dat.year.month$month <- as.integer(dat.year.month$month)
dat.year.month$year.month <- seq(nrow(dat.year.month))

# merge year.month table with population table to create year.month id
dat.inla <- merge(dat.inla,dat.year.month, by=c('year','month'))

# make sure that the order of the main data file matches that of the shapefile,
# otherwise the model will not be valid
dat.inla <- dat.inla[order(dat.inla$DRAWSEQ,dat.inla$sex,dat.inla$age,dat.inla$year.month),]

# add ID column for INLA
dat.inla$ID <- dat.inla$DRAWSEQ

# fix rownames
rownames(dat.inla) <- 1:nrow(dat.inla)

# variables for INLA model

dat.inla$year.month4 <- dat.inla$year.month3 <- dat.inla$year.month2 <- dat.inla$year.month
dat.inla$month4 <- dat.inla$month3 <- dat.inla$month2 <- dat.inla$month
dat.inla$ID3 <- dat.inla$ID2 <- dat.inla$ID
dat.inla$e <- 1:nrow(dat.inla)

# INLA

if(type==1){

# 1. Type I space-time interaction

fml1 <- deaths ~
	# global terms
        1 +                                                                     		# global intercept
        year.month +                                                            		# global slope
    	# month specific terms
        f(month, model="rw1",cyclic = TRUE) +                                   		# month specific intercept
        f(month2, year.month2, model="rw1",cyclic = TRUE) +                     		# month specific slope
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                      		# state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                        		# state specific slope (BYM)
	# random walk across time
        f(year.month3, model="rw1") +                                           		# rw1
        # overdispersion term
        f(e, model = "iid")                                                     		# overdispersion term

# INLA model

system.time(mod1 <-
    inla(formula = fml1,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type1_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod1$misc <- NULL
mod1$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod1,parameters.name)}
if(cluster==1){saveRDS(mod1,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type1_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod1 <- summary(mod1)
if(cluster==0){capture.output(inla.summary.mod1,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod1,file=paste0('../output/summary/',summary.name))}

# save plot of INLA model
#plot.name <- paste0('USA_rate_pred_type1_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod1 <- plot(mod1)
#plot(mod1, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod1$summary.fitted.values$mean,sd=mod1$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type1_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

if(type==2){

# 1. Type Ia space-time interaction

fml1a <- deaths ~
	# global terms
        1 +                                                                     		# global intercept
        year.month +                                                           			# global slope
    	# month specific terms
	f(month, model='rw1',cyclic = TRUE) +							# month specific intercept
	f(month2, year.month2, model='rw1', cyclic= TRUE) + 					# month specific slope
	# state-month specific terms
        f(month3, model="rw1",cyclic = TRUE,group=ID,control.group=list(model='besag',graph=USA.adj))+        		# state-month specific intercept (spatially-correlated)
        f(month4, year.month2, model="rw1",cyclic = TRUE,group=ID, control.group=list(model='besag',graph=USA.adj))+    # state-month specific slope (spatially-correlated)
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                      		# state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                        		# state specific slope (BYM)
	# random walk across time
        f(year.month3, model="rw1") +                                           		# rw1
        # overdispersion term
        f(e, model = "iid")                                                    		 	# overdispersion term

# INLA model

system.time(mod1a <-
    inla(formula = fml1a,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type1a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod1a$misc <- NULL
mod1a$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod1a,parameters.name)}
if(cluster==1){saveRDS(mod1a,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type1a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod1a <- summary(mod1a)
if(cluster==0){capture.output(inla.summary.mod1a,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod1a,file=paste0('../output/summary/',summary.name))}

# save plot of INLA model
#plot.name <- paste0('USA_rate_pred_type1a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod1a <- plot(mod1a)
#plot(mod1a, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod1a$summary.fitted.values$mean,sd=mod1a$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type1a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

if(type==3) {

# 2. Type II space-time interaction

fml2 <- deaths ~
	# global terms
        1 +                                                                     		# global intercept
        year.month +                                                            		# global slope
    	# month specific terms	
        f(month, model="rw1",cyclic = TRUE) +                                   		# month specific intercept
        f(month2, year.month2, model="rw1",cyclic = TRUE) +                     		# month specific slope
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                      		# state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                        		# state specific slope (BYM)
	# random walk across time
        f(year.month3, model="rw1") +                                           		# rw1
	# space-time interaction
        f(ID3,model="iid", group=year.month4,                                   		# type II space-time interaction
 	 control.group=list(model="rw1")) +
        #control.group=list(model="rw2")) +                                    
        #control.group=list(model="exchangeable")) + 
        # overdispersion term
        f(e, model = "iid")                                                     		# overdispersion term

# INLA model

system.time(mod2 <-
    inla(formula = fml2,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type2_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod2$misc <- NULL
mod2$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod2,parameters.name)}
if(cluster==1){saveRDS(mod2,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type2_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod2 <- summary(mod2)
if(cluster==0){capture.output(inla.summary.mod2,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod2,file=paste0('../output/summary/',summary.name))}

# save plot of INLA model
#plot.name <- paste0('USA_rate_pred_type2_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod2 <- plot(mod2)
#plot(mod2, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod2$summary.fitted.values$mean,sd=mod2$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type2_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

if(type==4) {

# 2. Type IIa space-time interaction

fml2a <- deaths ~
	# global terms
        1 +                                                                                     # global intercept
        year.month +                                                                            # global slope
     	# month specific terms
	f(month, model='rw1',cyclic = TRUE) +							# month specific intercept
	f(month2, year.month2, model='rw1', cyclic= TRUE) + 					# month specific slope
	# state-month specific terms
        f(month3, model="rw1",cyclic = TRUE,group=ID,control.group=list(model='besag',graph=USA.adj))+        		# state-month specific intercept (spatially-correlated)
        f(month4, year.month2, model="rw1",cyclic = TRUE,group=ID, control.group=list(model='besag',graph=USA.adj))+    # state-month specific slope (spatially-correlated)
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                                      # state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                                        # state specific slope (BYM)
	# random walk across time
        f(year.month3, model="rw1") +                                                           # rw1
	# space-time interaction
        f(ID3,model="iid", group=year.month4,                                                   # type II space-time interaction
 	 control.group=list(model="rw1")) +
        #control.group=list(model="rw2")) +                                    
        #control.group=list(model="exchangeable")) + 
        # overdispersion term
        f(e, model = "iid")                                                                     # overdispersion term

# INLA model

system.time(mod2a <-
    inla(formula = fml2a,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type2a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod2a$misc <- NULL
mod2a$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod2a,parameters.name)}
if(cluster==1){saveRDS(mod2a,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type2a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod2a <- summary(mod2a)
if(cluster==0){capture.output(inla.summary.mod2a,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod2a,file=paste0('../output/summary/',summary.name))}

# save plot of INLA model
#plot.name <- paste0('USA_rate_pred_type2a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod2a <- plot(mod2a)
#plot(mod2a, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod2a$summary.fitted.values$mean,sd=mod2a$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type2a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

if(type==5) {

# 3. Type III space-time interaction

fml3 <- deaths ~
	# global terms
        1 +                                                                     		# global intercept
        year.month +                                                           		 	# global slope
    	# month specific terms
        f(month, model="rw1",cyclic = TRUE) +                                   		# month specific intercept
        f(month2, year.month2, model="rw1",cyclic = TRUE) +                     		# month specific slope
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                      		# state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                        		# state specific slope (BYM)
	# random walk terms
        f(year.month3, model="rw1") +                                           		# rw1
	#f(year.month4,model="iid", group=ID3,				        		# type III space-time interaction
	f(year.month4,model="rw1", group=ID3,				        		# variation on model
	control.group=list(model="I")) +
	#control.group=list(model="besag", 							# variation on model  							
	#graph=USA.adj)) +									# variation on model 
        # overdispersion term
        f(e, model = "iid")                                                     		# overdispersion term

system.time(mod3 <-
    inla(formula = fml3,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    control.inla = list(strategy = "gaussian"),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type3_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod3$misc <- NULL
mod3$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod3,parameters.name)}
if(cluster==1){saveRDS(mod3,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type3_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod3 <- summary(mod3)
if(cluster==0){capture.output(inla.summary.mod3,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod3,file=paste0('../output/summary/',summary.name))}

# save plot of INLA model
#plot.name <- paste0('USA_rate_pred_type3_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod3 <- plot(mod3)
#plot(mod3, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod3$summary.fitted.values$mean,sd=mod3$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type3_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

if(type==6) {

# 3. Type IIIa space-time interaction

fml3 <- deaths ~
	# global terms
        1 +                                                                     		# global intercept
        year.month +                                                           		 	# global slope
     	# month specific terms
	f(month, model='rw1',cyclic = TRUE) +							# month specific intercept
	f(month2, year.month2, model='rw1', cyclic= TRUE) + 					# month specific slope
	# state-month specific terms
        f(month3, model="rw1",cyclic = TRUE,group=ID,control.group=list(model='besag',graph=USA.adj))+        		# state-month specific intercept (spatially-correlated)
        f(month4, year.month2, model="rw1",cyclic = TRUE,group=ID, control.group=list(model='besag',graph=USA.adj))+    # state-month specific slope (spatially-correlated)
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                      		# state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                        		# state specific slope (BYM)
	# random walks terms
        #f(year.month3, model="rw1") +                                           		# rw1 across time (put back?)
	#f(year.month4,model="iid", group=ID3,				        		# type III space-time interaction
	f(year.month4,model="rw1", group=ID3,				        		# variation on model
	control.group=list(model="I")) +
	#control.group=list(model="besag", 							# variation on model 							
	#graph=USA.adj)) +									# variation on model 
	# overdispersion term
        f(e, model = "iid")                                                     		# overdispersion term

system.time(mod3 <-
    inla(formula = fml3,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    control.inla = list(strategy = "gaussian"),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type3a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod3$misc <- NULL
mod3$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod3,parameters.name)}
if(cluster==1){saveRDS(mod3,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type3a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod3 <- summary(mod3)
if(cluster==0){capture.output(inla.summary.mod3,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod3,file=paste0('../output/summary/',summary.name))}

# save plot of INLA model
#plot.name <- paste0('USA_rate_pred_type3_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod3 <- plot(mod3)
#plot(mod3, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod3$summary.fitted.values$mean,sd=mod3$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type3a_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

if(type==7) {

# 4. Type IV space-time interaction

fml4 <- deaths ~
	# global terms
        1 +                                                                     		# global intercept
        year.month +                                                            		# global slope
    	# month specific terms
        f(month, model="rw1",cyclic = TRUE) +                                   		# month specific intercept
        f(month2, year.month2, model="rw1",cyclic = TRUE) +                     		# month specific slope
        # state specific terms
        f(ID, model="bym",graph=USA.adj) +                                     		 	# state specific intercept (BYM)
        f(ID2, year.month2, model="bym",graph=USA.adj) +                        		# state specific slope (BYM)
        # overdispersion term
        f(year.month3, model="rw1") +                                           		# rw1
        f(ID3,model="besag", graph=USA.adj,                                     		# type IV space-time interaction
        group=year.month4,
        control.group=list(model="rw1")) +
        f(e, model = "iid")                                                     		# overdispersion term

system.time(mod4 <-
    inla(formula = fml4,
    family = "poisson",
    data = dat.inla,
    E = pop.adj,
    control.compute = list(dic=TRUE),
    control.predictor = list(link = 1),
    control.inla = list(strategy = "gaussian"),
    verbose=TRUE
    ))

# save all parameters of INLA model
parameters.name <- paste0('USA_rate_pred_type4_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_parameters')
mod4$misc <- NULL
mod4$.args$.parent.frame <- NULL
if(cluster==0){saveRDS(mod4,parameters.name)}
if(cluster==1){saveRDS(mod4,paste0('../output/pred/',parameters.name))}

# save summary of INLA model
summary.name <- paste0('USA_rate_pred_type4_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.txt')
inla.summary.mod4 <- summary(mod4)
if(cluster==0){capture.output(inla.summary.mod4,file=summary.name)}
if(cluster==1){capture.output(inla.summary.mod4,file=paste0('../output/summary/',summary.name))}
# save plot of INLA model

#plot.name <- paste0('USA_rate_pred_type4_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end,'_summary.pdf')
#inla.plot.mod4 <- plot(mod4)
#plot(mod4, single=TRUE, pdf=TRUE, paper='a4r')

# save RDS of INLA results
plot.dat <- as.data.frame(cbind(dat.inla,rate.pred=mod4$summary.fitted.values$mean,sd=mod4$summary.fitted.values$sd))

# name of RDS output file then save
RDS.name <- paste0('USA_rate_pred_type4_',age,'_',sex.lookup[sex],'_',year.start,'_',year.end)
if(cluster==0){saveRDS(plot.dat,RDS.name)}
if(cluster==1){saveRDS(plot.dat,paste0('../output/pred/',RDS.name))}

}

sender <- "emailr349@gmail.com"
recipients <- c("r.parks15@imperial.ac.uk")
send.mail(from = sender,
          to = recipients,
          subject = paste0('Another age done'),
          body = "Well done for being you",
          smtp = list(host.name = "smtp.gmail.com", port = 465, 
                      user.name = "emailr349@gmail.com",            
                      passwd = "inlaisthebest", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

# this bracket ends the function at the top of the script
}

################

# input arguments into function to perform inference
#mapply(inla.function,age.sel=age.arg,sex.sel=sex.arg,year.start=year.start.arg,year.end=year.end.arg,type=type.arg,cluster=cluster.arg)

#mapply(inla.function,age.sel=c(75,85,0,5,15,25,35,45),sex.sel=1,year.start=1982,year.end=2010,type=2,cluster=0)
#mapply(inla.function,age.sel=c(55,65,75,85,0,5,15,25,35,45),sex.sel=2,year.start=1982,year.end=2010,type=2,cluster=0)

