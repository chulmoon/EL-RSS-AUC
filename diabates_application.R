library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)
library(NHANES) # diabates dataset

# parameter setting
# n = 20, 40, 80
# m = 2, 4, 5
params=expand.grid(nx=c(20,40,60,80),m=c(2,4,5))

#################################################################
############# functions for NHANES data application #############
#################################################################
main=function(pind,params,nsim=10000){
	source("./functions_BRSS.R")
	param=params[pind,]
	m=n=param$m
	k=l=param$nx/param$m
	nx=ny=param$nx
	
	len = cp = matrix(NA,nsim,3)
	rsse.rnew.ci = rssk.ci = srs.ci = matrix(NA,nsim,2)
	
	# load NHANES data
	nhanes.dat = NHANES
	
	# non-disease
	x = nhanes.dat %>%
		dplyr::filter(Diabetes=="No") %>%
		dplyr::select(BMI, DirectChol) %>%
		tidyr::drop_na()
	
	# disease
	y = nhanes.dat %>%
		dplyr::filter(Diabetes=="Yes") %>%
		dplyr::select(BMI, DirectChol) %>%
		tidyr::drop_na()
	
	xbmi = x$BMI
	ybmi = y$BMI
	
	# TRUE AUC
	AUC = 0
	for (i in 1:length(xbmi)){
		for (j in 1:length(ybmi)) {
			AUC=AUC+(ybmi[j]>=xbmi[i])
		}
	}
	AUC = AUC/(length(xbmi)*length(ybmi)) #0.733
	
	
	for (ii in 1:nsim) {
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			xx=x$BMI
			# concomitant of X
			cx=x$DirectChol
			
			# dist of Y
			yy=y$BMI
			# concomitant of y
			cy=y$DirectChol
			
			rssx=RSSampling::con.rss(xx,-cx,m=m,r=k)
			rssy=RSSampling::con.rss(yy,-cy,m=n,r=l)
			
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
############### diabetes data application ###################
#############################################################

# number of cores for parallel computation
numCores = parallel::detectCores()

cl=snow::makeCluster(numCores-1, type="SOCK")
diabetes.result=snow::parLapply(cl,1:nrow(params),main,params,nsim=5000)
snow::stopCluster(cl)
save(diabetes.result,file="./diabetes.result.Rdata")

