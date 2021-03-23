library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)

# parameter setup
# n_x = n_y = 20, 40, 80
# m = n = 2, 4, 5
# AUC = 0.6, 0.8, 0.9, 0.95
# judgement ranking  1, 0.9, 0.7

params=expand.grid(nx=c(20,40,80),m=c(2,4,5),
									 AUC=c(0.6,0.8,0.9,0.95),corxy=c(1,0.9,0.7))

#############################################################
################### Run BRSS simulations ####################
#############################################################
# number of cores for parallel computation
numCores = parallel::detectCores()

# normal
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_BRSS.R")
brss.normal.result=snow::parLapply(cl,1:nrow(params),
																	 main.brss.normal,params,nsim=5000)
snow::stopCluster(cl)
save(brss.normal.result,file="./brss.normal.result.Rdata")

# uniform
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_BRSS.R")
brss.uniform.result=snow::parLapply(cl,1:nrow(params),
																		main.brss.uniform,params,nsim=5000)
snow::stopCluster(cl)
save(brss.uniform.result,file="./brss.uniform.result.Rdata")

# lognormal
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_BRSS.R")
brss.lognormal.result=snow::parLapply(cl,1:nrow(params),
																			main.brss.lognormal,params,nsim=5000)
snow::stopCluster(cl)
save(brss.lognormal.result,file="./brss.lognormal.result.Rdata")
