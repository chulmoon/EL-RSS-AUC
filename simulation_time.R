library(emplik)
library(RSSampling)

# code for computation time

# parameters
# n = 20
# m = 2, 4, 5
# AUC = 0.6
# ranking 1
params=expand.grid(nx=c(20),m=c(2,4,5),
                   AUC=c(0.6),corxy=c(1))

########################################################
# SRS-EL
########################################################
main.srs=function(pind,params,nsim=10000){
  source("./functions_BRSS.R")
  param=params[pind,]
  m=n=param$m
  k=l=param$nx/param$m
  nx=ny=param$nx
  AUC=param$AUC
  corx=cory=param$corxy
  
  len = cp = matrix(NA,nsim,3)
  rsse.rnew.ci = rssk.ci = srs.ci = matrix(NA,nsim,2)
  
  for (sind in 1:nsim) {
    set.seed(10000*pind + sind)
    
    #############################################################
    # SRS
    srs.cond=TRUE
    while (srs.cond==T) {
      srsx=stats::rnorm(nx)
      srsy=stats::rnorm(ny,mu,sigma)
      
      srs.delhat = 0
      for (i in 1:length(srsx)){
        for (j in 1:length(srsy)) {
          srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
        }
      }
      srs.delhat = srs.delhat/(nx*ny)
      srs.cond=any(srs.delhat==1) # remove delta-hat=1 case
    }
    srs.el.res=srs.el(srsx,srsy,nx,ny)
    
    ## confidence interval
    srs.el.ci=c(srs.el.res$Low,srs.el.res$Up)
    ## converage probability
    srs.el.cp=(srs.el.ci[2]>=AUC)*(srs.el.ci[1]<=AUC)
    ## length
    srs.el.len=srs.el.res$Up-srs.el.res$Low
  }
}

# computation time for SRS-EL
ptm <- proc.time()
sim.res=main.srs(1,params,nsim=1000)
srs.el.time1=proc.time() - ptm

ptm <- proc.time()
sim.res=main.srs(2,params,nsim=1000)
srs.el.time2=proc.time() - ptm

ptm <- proc.time()
sim.res=main.srs(3,params,nsim=1000)
srs.el.time3=proc.time() - ptm


########################################################
# BRSS-EL
########################################################
main.brss.el=function(pind,params,nsim=10000){
  source("./functions_BRSS.R")
  param=params[pind,]
  m=n=param$m
  k=l=param$nx/param$m
  nx=ny=param$nx
  AUC=param$AUC
  corx=cory=param$corxy
  
  len = cp = matrix(NA,nsim,3)
  rsse.rnew.ci = rssk.ci = srs.ci = matrix(NA,nsim,2)
  
  for (sind in 1:nsim) {
    set.seed(10000*pind + sind)
    
    # RSS
    rss.cond=TRUE
    while(rss.cond==T){
      # dist of X
      nxsam = m^2*k
      x=stats::runif(nxsam)
      mux = 0.5
      sigma = sqrt(1/12)
      # concomitant of X
      cx=corx*((x-mux)/sigma) + sqrt(1-corx^2)*rnorm(nxsam)
      
      # dist of Y
      nysam= n^2*l
      theta = 1/(2-2*AUC)
      muy = theta/2
      sdy = sqrt(1/12*theta^2)
      y=stats::runif(nysam)*theta
      # concomitant of y
      cy=cory*((y-muy)/sdy) + sqrt(1-cory^2)*rnorm(nysam)
      
      rssx=con.rss.sim(x,cx,m=m,r=k)
      rssy=con.rss.sim(y,cy,m=n,r=l)
      
      rssx.vec=as.vector(rssx$sample.x)
      rssy.vec=as.vector(rssy$sample.x)
      
      delhat = 0
      for (i in 1:length(rssx.vec)){
        for (j in 1:length(rssy.vec)) {
          delhat=delhat+(rssy.vec[j]>=rssx.vec[i])
        }
      }
      delhat = delhat/(m*k*n*l)
      
      ######################################################
      # RSS + EL
      rss.el.rnew.res=rss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
      
      ## confidence interval
      rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
      ## coverage probability
      rss.el.rnew.cp=(rss.el.rnew.ci[2]>=AUC)*(rss.el.rnew.ci[1]<=AUC)
      ## length
      rss.el.rnew.len=rss.el.rnew.res$Up-rss.el.rnew.res$Low
      
      rss.cond=any( c(delhat==1, # remove delta-hat=1 case
                      is.na(rss.el.rnew.res$Low),
                      is.na(rss.el.rnew.res$Up) )
      )
    }
  }
}

# computation time for BRSS-EL
## set size m=n=2
ptm <- proc.time()
sim.res=main.brss.el(1,params,nsim=1000)
brss.el.time1=proc.time() - ptm

## set size m=n=4
ptm <- proc.time()
sim.res=main.brss.el(2,params,nsim=1000)
brss.el.time2=proc.time() - ptm

## set size m=n=5
ptm <- proc.time()
sim.res=main.brss.el(3,params,nsim=1000)
brss.el.time3=proc.time() - ptm

########################################################
# BRSS-KER
########################################################
main.brss.ker=function(pind,params,nsim=10000){
  source("./functions_BRSS.R")
  param=params[pind,]
  m=n=param$m
  k=l=param$nx/param$m
  nx=ny=param$nx
  AUC=param$AUC
  corx=cory=param$corxy
  
  len = cp = matrix(NA,nsim,3)
  rsse.rnew.ci = rssk.ci = srs.ci = matrix(NA,nsim,2)
  
  for (sind in 1:nsim) {
    set.seed(10000*pind + sind)
    
    # RSS
    rss.cond=TRUE
    while(rss.cond==T){
      # dist of X
      nxsam = m^2*k
      x=stats::runif(nxsam)
      mux = 0.5
      sigma = sqrt(1/12)
      # concomitant of X
      cx=corx*((x-mux)/sigma) + sqrt(1-corx^2)*rnorm(nxsam)
      
      # dist of Y
      nysam= n^2*l
      theta = 1/(2-2*AUC)
      muy = theta/2
      sdy = sqrt(1/12*theta^2)
      y=stats::runif(nysam)*theta
      # concomitant of y
      cy=cory*((y-muy)/sdy) + sqrt(1-cory^2)*rnorm(nysam)
      
      rssx=con.rss.sim(x,cx,m=m,r=k)
      rssy=con.rss.sim(y,cy,m=n,r=l)
      
      rssx.vec=as.vector(rssx$sample.x)
      rssy.vec=as.vector(rssy$sample.x)
      
      delhat = 0
      for (i in 1:length(rssx.vec)){
        for (j in 1:length(rssy.vec)) {
          delhat=delhat+(rssy.vec[j]>=rssx.vec[i])
        }
      }
      delhat = delhat/(m*k*n*l)
      
      ######################################################
      # RSS + Kernel
      rss.ker.res=rss.ker(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
      
      ## confidence interval
      rss.ker.ci=rss.ker.res
      ## coverage probability
      rss.ker.cp=(rss.ker.ci[2]>=AUC)*(rss.ker.ci[1]<=AUC)
      ## length
      rss.ker.len=rss.ker.res[2]-rss.ker.res[1]
      
      rss.cond=any( c(delhat==1, # remove delta-hat=1 case
                      is.na(rss.ker.ci[1]),
                      is.na(rss.ker.ci[2]))
      )
      
    }
  }
}

# computation time for BRSS-KER
## set size m=n=2
ptm <- proc.time()
sim.res=main.brss.ker(1,params,nsim=1000)
brss.ker.time1=proc.time() - ptm

## set size m=n=4
ptm <- proc.time()
sim.res=main.brss.ker(2,params,nsim=1000)
brss.ker.time2=proc.time() - ptm

## set size m=n=5
ptm <- proc.time()
sim.res=main.brss.ker(3,params,nsim=1000)
brss.ker.time3=proc.time() - ptm

