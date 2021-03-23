
#############################################################
############## main URSS simulation functions ###############
#############################################################

########################################################
# normal
########################################################
main.urss.normal=function(pind,params,nsim=10000){
	source("./functions_URSS.R")
	param=params[pind,]
	m=n=param$m
	nx=ny=param$nx
	prop=param$prop
	ki=c(nx*0.5,nx*0.5) # number of cycles of X of URSS
	lr=c(ny*prop,ny*(1-prop)) # number of cycles of Y of URSS
	
	AUC=param$AUC
	corx=cory=param$corxy
	
	len = cp = matrix(NA,nsim,1)
	rsse.rnew.ci = matrix(NA,nsim,2)
	
	for (sind in 1:nsim) {
		set.seed(70000*pind + sind)
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			# dist of X
			nxsam = 10000
			x=rnorm(nxsam)
			# concomitant of X
			cx=corx*x + sqrt(1-corx^2)*rnorm(nxsam)
			
			# dist of Y
			ny.sam= 10000
			mu = sqrt(5)*qnorm(AUC)
			sigma = 2
			y=rnorm(ny.sam,mu,sigma)
			# concomitant of y
			cy=cory*((y-mu)/sigma) + sqrt(1-cory^2)*rnorm(ny.sam)
			
			rssx=list()
			rssx.temp=RSSampling::con.rss(x,cx,m=m,r=max(ki))
			for (kk in 1:m){
				rssx[[kk]]=rssx.temp$sample.x[1:ki[kk],kk]
			}
			rssy=list()
			rssy.temp=RSSampling::con.rss(y,cy,m=n,r=max(lr))
			for (ll in 1:n){
				rssy[[ll]]=rssy.temp$sample.x[1:lr[ll],ll]
			}
			
			rssx.vec=unlist(rssx)
			rssy.vec=unlist(rssy)
			
			## delta hat
			delhat = 0
			delhatx = rep(NA,length(rssx))
			delhaty = rep(NA,length(rssy))
			for (ii in 1:length(rssx)){
				delhatxij = rep(NA,length(rssx[[ii]]))
				for (jj in 1:length(rssx[[ii]])){
					for (rr in 1:length(rssy)) {
						delhaty[rr]=mean(rssy[[rr]]>=rssx[[ii]][jj])
					}
					delhatxij[jj]=mean(delhaty)
				}
				delhatx[ii]=mean(delhatxij)
			}
			delhat=mean(delhatx)
			
			rss.el.rnew.res=urss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,ki,n,lr,nx,ny)
			rss.cond=any( c(delhat==1,is.na(rss.el.rnew.res$Low),is.na(rss.el.rnew.res$Up)) )
		}
		
		#############################################################
		## confidence interval
		rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
		## coverage probability
		rss.el.rnew.cp=(rss.el.rnew.ci[2]>=AUC)*(rss.el.rnew.ci[1]<=AUC)
		## length
		rss.el.rnew.len=rss.el.rnew.res$Up-rss.el.rnew.res$Low
		
		# save results
		## confidence interval
		rsse.rnew.ci[sind,]=rss.el.rnew.ci
		## coverage probability
		cp[sind,]=c(rss.el.rnew.cp)
		## length
		len[sind,]=c(rss.el.rnew.len)
	}
	# return results
	ciavg=apply(rsse.rnew.ci,2,mean)
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}

########################################################
# uniform
########################################################
main.urss.uniform=function(pind,params,nsim=10000){
	source("./functions_URSS.R")
	param=params[pind,]
	m=n=param$m
	nx=ny=param$nx
	prop=param$prop
	ki=c(nx*0.5,nx*0.5) # number of cycles of X of URSS
	lr=c(ny*prop,ny*(1-prop)) # number of cycles of Y of URSS
	
	AUC=param$AUC
	corx=cory=param$corxy
	
	len = cp = matrix(NA,nsim,1)
	rsse.rnew.ci = matrix(NA,nsim,2)
	
	for (sind in 1:nsim) {
		set.seed(110000*pind + sind)
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			# dist of X
			nxsam = 10000
			x=runif(nxsam)
			mux = 0.5
			sigma = sqrt(1/12)
			# concomitant of X
			cx=corx*((x-mux)/sigma) + sqrt(1-corx^2)*rnorm(nxsam)
			
			# dist of Y
			ny.sam= 10000
			theta = 1/(2-2*AUC)
			muy = theta/2
			sdy = sqrt(1/12*theta^2)
			y=runif(ny.sam)*theta
			# concomitant of y
			cy=cory*((y-muy)/sdy) + sqrt(1-cory^2)*rnorm(ny.sam)
			
			rssx=list()
			rssx.temp=RSSampling::con.rss(x,cx,m=m,r=max(ki))
			for (kk in 1:m){
				rssx[[kk]]=rssx.temp$sample.x[1:ki[kk],kk]
			}
			rssy=list()
			rssy.temp=RSSampling::con.rss(y,cy,m=n,r=max(lr))
			for (ll in 1:n){
				rssy[[ll]]=rssy.temp$sample.x[1:lr[ll],ll]
			}
			
			rssx.vec=unlist(rssx)
			rssy.vec=unlist(rssy)
			
			## delta hat
			delhat = 0
			delhatx = rep(NA,length(rssx))
			delhaty = rep(NA,length(rssy))
			for (ii in 1:length(rssx)){
				delhatxij = rep(NA,length(rssx[[ii]]))
				for (jj in 1:length(rssx[[ii]])){
					for (rr in 1:length(rssy)) {
						delhaty[rr]=mean(rssy[[rr]]>=rssx[[ii]][jj])
					}
					delhatxij[jj]=mean(delhaty)
				}
				delhatx[ii]=mean(delhatxij)
			}
			delhat=mean(delhatx)
			
			rss.el.rnew.res=urss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,ki,n,lr,nx,ny)
			rss.cond=any( c(delhat==1,is.na(rss.el.rnew.res$Low),is.na(rss.el.rnew.res$Up)) )
		}
		
		#############################################################
		## confidence interval
		rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
		## coverage probability
		rss.el.rnew.cp=(rss.el.rnew.ci[2]>=AUC)*(rss.el.rnew.ci[1]<=AUC)
		## length
		rss.el.rnew.len=rss.el.rnew.res$Up-rss.el.rnew.res$Low
		
		# save results
		## confidence interval
		rsse.rnew.ci[sind,]=rss.el.rnew.ci
		## coverage probability
		cp[sind,]=c(rss.el.rnew.cp)
		## length
		len[sind,]=c(rss.el.rnew.len)
	}
	# return results
	ciavg=apply(rsse.rnew.ci,2,mean)
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}


########################################################
# lognormal
########################################################
main.urss.lognormal=function(pind,params,nsim=10000){
	source("./functions_URSS.R")
	param=params[pind,]
	m=n=param$m
	nx=ny=param$nx
	prop=param$prop
	ki=c(nx*0.5,nx*0.5) # number of cycles of X of URSS
	lr=c(ny*prop,ny*(1-prop)) # number of cycles of Y of URSS
	
	AUC=param$AUC
	corx=cory=param$corxy
	
	len = cp = matrix(NA,nsim,1)
	rsse.rnew.ci = matrix(NA,nsim,2)
	
	for (sind in 1:nsim) {
		set.seed(90000*pind + sind)
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			# dist of X
			nxsam = 10000
			x=rlnorm(nxsam)
			# concomitant of X
			mux=exp(1/2)
			sdx=sqrt( (exp(1)-1)*exp(1) )
			cx=corx*((x-mux)/sdx) + sqrt(1-corx^2)*rnorm(nxsam)
			
			# dist of Y
			ny.sam= 10000
			sig=2
			muy.temp = sqrt(sig^2+1)*qnorm(AUC)
			sdy.temp = sig
			muy=exp(muy.temp+(sdy.temp^2)/2)
			sdy=sqrt( (exp(sdy.temp^2)-1)*exp(2*muy.temp+sdy.temp^2) )
			y=rlnorm(ny.sam,meanlog=muy.temp,sdlog=sdy.temp)
			# concomitant of y
			cy=cory*((y-muy)/sdy) + sqrt(1-cory^2)*rnorm(ny.sam)
			
			rssx=list()
			rssx.temp=RSSampling::con.rss(x,cx,m=m,r=max(ki))
			for (kk in 1:m){
				rssx[[kk]]=rssx.temp$sample.x[1:ki[kk],kk]
			}
			rssy=list()
			rssy.temp=RSSampling::con.rss(y,cy,m=n,r=max(lr))
			for (ll in 1:n){
				rssy[[ll]]=rssy.temp$sample.x[1:lr[ll],ll]
			}
			
			rssx.vec=unlist(rssx)
			rssy.vec=unlist(rssy)
			
			## delta hat
			delhat = 0
			delhatx = rep(NA,length(rssx))
			delhaty = rep(NA,length(rssy))
			for (ii in 1:length(rssx)){
				delhatxij = rep(NA,length(rssx[[ii]]))
				for (jj in 1:length(rssx[[ii]])){
					for (rr in 1:length(rssy)) {
						delhaty[rr]=mean(rssy[[rr]]>=rssx[[ii]][jj])
					}
					delhatxij[jj]=mean(delhaty)
				}
				delhatx[ii]=mean(delhatxij)
			}
			delhat=mean(delhatx)
			
			rss.el.rnew.res=urss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,ki,n,lr,nx,ny)
			rss.cond=any( c(delhat==1,is.na(rss.el.rnew.res$Low),is.na(rss.el.rnew.res$Up)) )
		}
		
		#############################################################
		## confidence interval
		rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
		## coverage probability
		rss.el.rnew.cp=(rss.el.rnew.ci[2]>=AUC)*(rss.el.rnew.ci[1]<=AUC)
		## length
		rss.el.rnew.len=rss.el.rnew.res$Up-rss.el.rnew.res$Low
		
		# save results
		## confidence interval
		rsse.rnew.ci[sind,]=rss.el.rnew.ci
		## coverage probability
		cp[sind,]=c(rss.el.rnew.cp)
		## length
		len[sind,]=c(rss.el.rnew.len)
	}
	# return results
	ciavg=apply(rsse.rnew.ci,2,mean)
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}

#############################################################
############## main SRS simulation functions ###############
#############################################################

########################################################
# normal
########################################################
main.srs.normal=function(jj,params,nsim=10000){
	param=params[jj,]
	m=n=param$m
	nx=ny=param$nx
	prop=param$prop
	ki=c(nx*prop,nx*(1-prop))
	lr=c(ny*prop,ny*(1-prop))
	
	AUC=param$AUC
	corx=cory=param$corxy
	
	len = cp = matrix(NA,nsim,1)
	srs.ci = matrix(NA,nsim,2)
	
	# adjusted EL
	el.auc <- function(mu,x,r.adj) {
		elres = emplik::el.test(x,mu)
		# adjust log EL
		elres$"-2LLR" = elres$"-2LLR"*r.adj
		return(elres)
	}
	
	#####################################
	# SRS + EL
	srsv10 = function(xij,y,ny) {
		return(1/(ny)*sum(xij<=y))
	}
	srsv01 = function(yrs,x,nx) {
		return(1/(nx)*sum(yrs>=x))
	}
	
	srs.el = function(srsx,srsy,nx,ny){
		## One minus U
		srsOU = rep(NA,length(srsy))
		for (j in 1:length(srsy)) {
			srsOU[j] = (sum(srsy[j]>=srsx)/(length(srsx)))
		}
		
		## delta hat
		srs.delhat = 0
		for (i in 1:length(srsx)){
			for (j in 1:length(srsy)) {
				srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
			}
		}
		srs.delhat = srs.delhat/(nx*ny)
		
		## sum(1-Uhat-delhat)^2
		srsOUDsq = 0
		for (j in 1:length(srsy)) {
			srsOUDsq = srsOUDsq + (sum(srsy[j]>=srsx)/(length(srsx))-srs.delhat)^2
		}
		
		## S^2
		srsS10sq=1/(nx-1)*( sum((purrr::map_dbl(srsx,srsv10,srsy,nx)-srs.delhat)^2) )
		srsS01sq=1/(ny-1)*( sum((purrr::map_dbl(srsy,srsv01,srsx,ny)-srs.delhat)^2) )
		srsSsq = (nx*srsS01sq+ny*srsS10sq)/(nx+ny)
		
		## adjustment
		srs.r.adj = (nx)/(nx+ny)*srsOUDsq/(ny*srsSsq)
		srs.el.res=emplik::findUL2(step=0.005, fun=el.auc, MLE=srs.delhat, x=srsOU, r.adj=srs.r.adj)
	}
	
	
	for (ii in 1:nsim) {
		set.seed(70000*jj + ii)
		
		mu = sqrt(5)*qnorm(AUC)
		sigma = 2
		
		# SRS
		srs.cond=TRUE
		while (srs.cond==T) {
			srsx=rnorm(nx)
			srsy=rnorm(ny,mu,sigma)
			
			srs.delhat = 0
			for (i in 1:length(srsx)){
				for (j in 1:length(srsy)) {
					srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
				}
			}
			srs.delhat = srs.delhat/(nx*ny)
			srs.cond=any(c(srs.delhat==1))
		}
		
		######################################################
		# SRS + EL
		srs.el.res=srs.el(srsx,srsy,nx,ny)
		srs.el.ci=c(srs.el.res$Low,srs.el.res$Up)
		srs.el.cp=(srs.el.ci[2]>=AUC)*(srs.el.ci[1]<=AUC)
		srs.el.len=srs.el.res$Up-srs.el.res$Low
		
		# results
		srs.ci[ii,]=srs.el.ci
		cp[ii,]=srs.el.cp
		len[ii,]=srs.el.len
	}
	
	ciavg=rbind(apply(srs.ci,2,mean))
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}


########################################################
# uniform
########################################################
main.srs.uniform=function(jj,params,nsim=10000){
	param=params[jj,]
	m=n=param$m
	nx=ny=param$nx
	prop=param$prop
	ki=c(nx*prop,nx*(1-prop))
	lr=c(ny*prop,ny*(1-prop))
	
	AUC=param$AUC
	corx=cory=param$corxy
	
	len = cp = matrix(NA,nsim,1)
	srs.ci = matrix(NA,nsim,2)
	
	# adjusted EL
	el.auc <- function(mu,x,r.adj) {
		elres = emplik::el.test(x,mu)
		# adjust log EL
		elres$"-2LLR" = elres$"-2LLR"*r.adj
		return(elres)
	}
	
	
	#####################################
	# SRS + EL
	srsv10 = function(xij,y,ny) {
		return(1/(ny)*sum(xij<=y))
	}
	srsv01 = function(yrs,x,nx) {
		return(1/(nx)*sum(yrs>=x))
	}
	
	srs.el = function(srsx,srsy,nx,ny){
		## One minus U
		srsOU = rep(NA,length(srsy))
		for (j in 1:length(srsy)) {
			srsOU[j] = (sum(srsy[j]>=srsx)/(length(srsx)))
		}
		
		## delta hat
		srs.delhat = 0
		for (i in 1:length(srsx)){
			for (j in 1:length(srsy)) {
				srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
			}
		}
		srs.delhat = srs.delhat/(nx*ny)
		
		## sum(1-Uhat-delhat)^2
		srsOUDsq = 0
		for (j in 1:length(srsy)) {
			srsOUDsq = srsOUDsq + (sum(srsy[j]>=srsx)/(length(srsx))-srs.delhat)^2
		}
		
		## S^2
		srsS10sq=1/(nx-1)*( sum((purrr::map_dbl(srsx,srsv10,srsy,nx)-srs.delhat)^2) )
		srsS01sq=1/(ny-1)*( sum((purrr::map_dbl(srsy,srsv01,srsx,ny)-srs.delhat)^2) )
		srsSsq = (nx*srsS01sq+ny*srsS10sq)/(nx+ny)
		
		## adjustment
		srs.r.adj = (nx)/(nx+ny)*srsOUDsq/(ny*srsSsq)
		srs.el.res=emplik::findUL2(step=0.005, fun=el.auc, MLE=srs.delhat, x=srsOU, r.adj=srs.r.adj)
	}
	
	
	for (ii in 1:nsim) {
		set.seed(110000*jj + ii)
		
		theta = 1/(2-2*AUC)
		
		# SRS
		srs.cond=TRUE
		while (srs.cond==T) {
			srsx=runif(nx)
			srsy=runif(ny)*theta
			
			srs.delhat = 0
			for (i in 1:length(srsx)){
				for (j in 1:length(srsy)) {
					srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
				}
			}
			srs.delhat = srs.delhat/(nx*ny)
			srs.cond=any(c(srs.delhat==1))
		}
		
		######################################################
		# SRS + EL
		srs.el.res=srs.el(srsx,srsy,nx,ny)
		srs.el.ci=c(srs.el.res$Low,srs.el.res$Up)
		srs.el.cp=(srs.el.ci[2]>=AUC)*(srs.el.ci[1]<=AUC)
		srs.el.len=srs.el.res$Up-srs.el.res$Low
		
		# results
		srs.ci[ii,]=srs.el.ci
		cp[ii,]=srs.el.cp
		len[ii,]=srs.el.len
	}
	ciavg=rbind(apply(srs.ci,2,mean))
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}


########################################################
# lognormal
########################################################
main.srs.lognormal=function(jj,params,nsim=10000){
	param=params[jj,]
	m=n=param$m
	nx=ny=param$nx
	prop=param$prop
	ki=c(nx*prop,nx*(1-prop))
	lr=c(ny*prop,ny*(1-prop))
	
	AUC=param$AUC
	corx=cory=param$corxy
	
	len = cp = matrix(NA,nsim,1)
	srs.ci = matrix(NA,nsim,2)
	
	# adjusted EL
	el.auc <- function(mu,x,r.adj) {
		elres = emplik::el.test(x,mu)
		# adjust log EL
		elres$"-2LLR" = elres$"-2LLR"*r.adj
		return(elres)
	}
	
	#####################################
	# SRS + EL
	srsv10 = function(xij,y,ny) {
		return(1/(ny)*sum(xij<=y))
	}
	srsv01 = function(yrs,x,nx) {
		return(1/(nx)*sum(yrs>=x))
	}
	
	srs.el = function(srsx,srsy,nx,ny){
		## One minus U
		srsOU = rep(NA,length(srsy))
		for (j in 1:length(srsy)) {
			srsOU[j] = (sum(srsy[j]>=srsx)/(length(srsx)))
		}
		
		## delta hat
		srs.delhat = 0
		for (i in 1:length(srsx)){
			for (j in 1:length(srsy)) {
				srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
			}
		}
		srs.delhat = srs.delhat/(nx*ny)
		
		## sum(1-Uhat-delhat)^2
		srsOUDsq = 0
		for (j in 1:length(srsy)) {
			srsOUDsq = srsOUDsq + (sum(srsy[j]>=srsx)/(length(srsx))-srs.delhat)^2
		}
		
		## S^2
		srsS10sq=1/(nx-1)*( sum((purrr::map_dbl(srsx,srsv10,srsy,nx)-srs.delhat)^2) )
		srsS01sq=1/(ny-1)*( sum((purrr::map_dbl(srsy,srsv01,srsx,ny)-srs.delhat)^2) )
		srsSsq = (nx*srsS01sq+ny*srsS10sq)/(nx+ny)
		
		## adjustment
		srs.r.adj = (nx)/(nx+ny)*srsOUDsq/(ny*srsSsq)
		srs.el.res=emplik::findUL2(step=0.005, fun=el.auc, MLE=srs.delhat, x=srsOU, r.adj=srs.r.adj)
	}
	
	
	for (ii in 1:nsim) {
		set.seed(90000*jj + ii)
		
		muy.temp = sqrt(101)*qnorm(AUC)
		sdy.temp = 10
		
		# SRS
		srs.cond=TRUE
		while (srs.cond==T) {
			srsx=rlnorm(nx)
			srsy=rlnorm(ny,meanlog=muy.temp,sdlog=sdy.temp)
			
			srs.delhat = 0
			for (i in 1:length(srsx)){
				for (j in 1:length(srsy)) {
					srs.delhat=srs.delhat+(srsy[j]>=srsx[i])
				}
			}
			srs.delhat = srs.delhat/(nx*ny)
			srs.cond=any(c(srs.delhat==1))
		}
		
		######################################################
		# SRS + EL
		srs.el.res=srs.el(srsx,srsy,nx,ny)
		srs.el.ci=c(srs.el.res$Low,srs.el.res$Up)
		srs.el.cp=(srs.el.ci[2]>=AUC)*(srs.el.ci[1]<=AUC)
		srs.el.len=srs.el.res$Up-srs.el.res$Low
		
		# results
		srs.ci[ii,]=srs.el.ci
		cp[ii,]=srs.el.cp
		len[ii,]=srs.el.len
	}
	ciavg=rbind(apply(srs.ci,2,mean))
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}
