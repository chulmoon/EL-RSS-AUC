library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)

# parameter setup
# n_x = {40, 100, 200}
# n_y = 40
# m = n = 2
# AUC = 0.8
# judgment ranking  1

params=expand.grid(nx=c(40,100,200),ny=40,m=2,n=2,
                   AUC=c(0.8),corxy=c(1))

#############################################################
################### Run BRSS simulations ####################
#############################################################

# number of cores for parallel computation
numCores = parallel::detectCores()

cl=snow::makeCluster(numCores-1, type="SOCK")
source("./simulation_functions_BRSS.R")
brss.normal.result=snow::parLapply(cl,1:nrow(params),
                                   main.brss.normal.reference,params,nsim=1000)
snow::stopCluster(cl)
save(brss.normal.result,file="./brss.normal.reference.result.Rdata")
