library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)

# parameter setup
## n_x = n_y = 20, 40, 80
## m = n = 2
## p_y = 0.3, 0.4, 0.5, 0.6, 0.7
## AUC = 0.6, 0.8, 0.9, 0.95
## judgement ranking rho_x = rho_y = 1

params=expand.grid(nx=c(20,40,80),m=c(2),prop=c(0.3,0.4,0.5,0.6,0.7),
									 AUC=c(0.6,0.8,0.9,0.95),corxy=c(1))

#############################################################
################### Run URSS simulations ####################
#############################################################
# number of cores for parallel computation
numCores = parallel::detectCores()

# normal
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_URSS.R")
urss.normal.result=snow::parLapply(cl,1:nrow(params),
																	 main.urss.normal,params,nsim=5000)
snow::stopCluster(cl)
save(urss.normal.result,file="./urss.normal.result.Rdata")

# uniform
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_URSS.R")
urss.uniform.result=snow::parLapply(cl,1:nrow(params),
																		main.urss.uniform,params,nsim=5000)
snow::stopCluster(cl)
save(urss.uniform.result,file="./urss.uniform.result.Rdata")

# lognormal
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_URSS.R")
urss.lognormal.result=snow::parLapply(cl,1:nrow(params),
																			main.urss.lognormal,params,nsim=5000)
snow::stopCluster(cl)
save(urss.lognormal.result,file="./urss.lognormal.result.Rdata")


#############################################################
################### Run SRS simulations #####################
#############################################################
# normal
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_URSS.R")
urss.srs.normal.result=snow::parLapply(cl,1:nrow(params),
																			 main.srs.normal,params,nsim=5000)
snow::stopCluster(cl)
save(urss.srs.normal.result,file="./urss.srs.normal.result.Rdata")

# uniform
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_URSS.R")
urss.srs.uniform.result=snow::parLapply(cl,1:nrow(params),
																				main.srs.uniform,params,nsim=5000)
snow::stopCluster(cl)
save(urss.srs.uniform.result,file="./urss.srs.uniform.result.Rdata")

# lognormal
cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_URSS.R")
urss.srs.lognormal.result=snow::parLapply(cl,1:nrow(params),
																					main.srs.lognormal,params,nsim=5000)
snow::stopCluster(cl)
save(urss.srs.lognormal.result,file="./urss.srs.lognormal.result.Rdata")