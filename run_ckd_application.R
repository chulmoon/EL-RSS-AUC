library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)

# parameter setting
# n = 20, 40, 60, 80
# m = 2, 4
params=expand.grid(nx=c(20,40,60,80),m=c(2,4))

#################################################################
############## functions for CKD data application ###############
#################################################################

main=function(pind,params,nsim=10000){
  source("./functions_BRSS.R")
  param=params[pind,]
  m=n=param$m
  k=l=param$nx/param$m
  nx=ny=param$nx
  
  len = cp = matrix(NA,nsim,3)
  rsse.rnew.ci = rssk.ci = srs.ci = matrix(NA,nsim,2)
  
  # load data
  load("./data_ckd.Rdata")
  # non-disease
  x = data.ckd %>%
    dplyr::filter(CKD==0)
  
  # disease
  y = data.ckd %>%
    dplyr::filter(CKD==1)
  
  xx=-x$GFR
  # concomitant of X
  cx=x$RIDAGEYR
  
  # dist of Y
  yy=-y$GFR
  # concomitant of y
  cy=y$RIDAGEYR
  
  
  # TRUE AUC
  AUC = 0
  for (i in 1:length(xx)){
    for (j in 1:length(yy)) {
      AUC=AUC+(yy[j]>=xx[i])
    }
  }
  AUC = AUC/(length(xx)*length(yy))
  
  
  for (ii in 1:nsim) {
    set.seed(70000*pind + ii)
    # RSS
    rss.cond=TRUE
    while(rss.cond==T){
      rssx=RSSampling::con.rss(xx,cx,m=m,r=k)
      rssy=RSSampling::con.rss(yy,cy,m=n,r=l)
      
      rssx.vec=as.vector(rssx$sample.x)
      rssy.vec=as.vector(rssy$sample.x)
      
      delhat = 0
      for (i in 1:length(rssx.vec)){
        for (j in 1:length(rssy.vec)) {
          delhat=delhat+(rssy.vec[j]>=rssx.vec[i])
        }
      }
      delhat = delhat/(m*k*n*l)
      rss.el.rnew.res=rss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
      rss.cond=any( c(delhat==1,is.na(rss.el.rnew.res$Low),is.na(rss.el.rnew.res$Up)) )
    }
    
    
    # SRS
    srs.cond=TRUE
    while (srs.cond==T) {
      srsx=sample(xx,nx)
      srsy=sample(yy,ny)
      
      srs.delhat = 0
      for (i in 1:length(srsx)){
        for (j in 1:length(srsy)) {
          srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
        }
      }
      srs.delhat = srs.delhat/(nx*ny)
      srs.cond=any(c(srs.delhat==1))
    }
    
    #############################################################
    # RSS + EL
    rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
    rss.el.rnew.cp=(rss.el.rnew.ci[2]>=AUC)*(rss.el.rnew.ci[1]<=AUC)
    rss.el.rnew.len=rss.el.rnew.res$Up-rss.el.rnew.res$Low
    
    ######################################################
    # RSS + Kernel
    
    rss.ker.res=rss.ker(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
    rss.ker.ci=rss.ker.res
    rss.ker.cp=(rss.ker.ci[2]>=AUC)*(rss.ker.ci[1]<=AUC)
    rss.ker.len=rss.ker.res[2]-rss.ker.res[1]
    
    ######################################################
    # SRS + EL
    srs.el.res=srs.el(srsx,srsy,nx,ny)
    srs.el.ci=c(srs.el.res$Low,srs.el.res$Up)
    srs.el.cp=(srs.el.ci[2]>=AUC)*(srs.el.ci[1]<=AUC)
    srs.el.len=srs.el.res$Up-srs.el.res$Low
    
    # results
    ## confidence interval
    rsse.rnew.ci[ii,]=rss.el.rnew.ci
    rssk.ci[ii,]=rss.ker.ci
    srs.ci[ii,]=srs.el.ci
    ## coverage probability
    cp[ii,]=c(rss.el.rnew.cp, rss.ker.cp, srs.el.cp)
    ## length
    len[ii,]=c(rss.el.rnew.len, rss.ker.len, srs.el.len)
  }
  ciavg=rbind(
    apply(rsse.rnew.ci,2,mean),
    apply(rssk.ci,2,mean),
    apply(srs.ci,2,mean)
  )
  cpavg=apply(cp,2,mean)
  lenavg=apply(len,2,mean)
  return(list(ci=ciavg,cp=cpavg,len=lenavg))
}


#############################################################
################### CKD data application ####################
#############################################################

# number of cores for parallel computation
numCores = parallel::detectCores()

cl=snow::makeCluster(numCores-1, type="SOCK")
sim.res=snow::parLapply(cl,1:nrow(params),main,params,nsim=5000)
snow::stopCluster(cl)

save(sim.res,file="./sim.ckd.res.Rdata")
