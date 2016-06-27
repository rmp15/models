rm(list=ls())

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
month.lookup <- c('January','February','March','April','May','June','July',
		'August','September','October','November','December')

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



