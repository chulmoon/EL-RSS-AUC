
#############################################################
############## main BRSS simulation functions ###############
#############################################################

########################################################
# normal
########################################################
main.brss.normal=function(pind,params,nsim=10000){
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
		print(sind)
		#############################################################
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			# dist of X
			nxsam = 10000
			x=stats::rnorm(nxsam)
			# concomitant of X
			cx=corx*x + sqrt(1-corx^2)*rnorm(nxsam)
			
			# dist of Y
			ny.sam= 10000
			mu = sqrt(5)*qnorm(AUC)
			sigma = 2
			y=stats::rnorm(ny.sam,mu,sigma)
			# concomitant of y
			cy=cory*((y-mu)/sigma) + sqrt(1-cory^2)*rnorm(ny.sam)
			
			rssx=RSSampling::con.rss(x,cx,m=m,r=k)
			rssy=RSSampling::con.rss(y,cy,m=n,r=l)
			
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
											is.na(rss.el.rnew.res$Low),
											is.na(rss.el.rnew.res$Up),
										  is.na(rss.ker.ci[1]),
										  is.na(rss.ker.ci[2]))
										)
		}
		
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
		
		# save results
		## confidence interval
		rsse.rnew.ci[sind,]=rss.el.rnew.ci
		rssk.ci[sind,]=rss.ker.ci
		srs.ci[sind,]=srs.el.ci
		## coverage probability
		cp[sind,]=c(rss.el.rnew.cp, rss.ker.cp, srs.el.cp)
		len[sind,]=c(rss.el.rnew.len, rss.ker.len, srs.el.len)
	}
	# return results
	ciavg=rbind(
		apply(rsse.rnew.ci,2,mean),
		apply(rssk.ci,2,mean),
		apply(srs.ci,2,mean)
	)
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}

########################################################
# uniform
########################################################
main.brss.uniform=function(pind,params,nsim=10000){
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
		set.seed(20000*pind + sind)
		#############################################################
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			# dist of X
			nxsam = 10000
			x=stats::runif(nxsam)
			mux = 0.5
			sigma = sqrt(1/12)
			# concomitant of X
			cx=corx*((x-mux)/sigma) + sqrt(1-corx^2)*rnorm(nxsam)
			
			# dist of Y
			ny.sam= 10000
			theta = 1/(2-2*AUC)
			muy = theta/2
			sdy = sqrt(1/12*theta^2)
			y=stats::runif(ny.sam)*theta
			# concomitant of y
			cy=cory*((y-muy)/sdy) + sqrt(1-cory^2)*rnorm(ny.sam)
			
			rssx=RSSampling::con.rss(x,cx,m=m,r=k)
			rssy=RSSampling::con.rss(y,cy,m=n,r=l)
			
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
											is.na(rss.el.rnew.res$Low),
											is.na(rss.el.rnew.res$Up),
											is.na(rss.ker.ci[1]),
											is.na(rss.ker.ci[2]))
			)
		}
		
		#############################################################
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
			srs.cond=any(srs.delhat==1) # remove delta-hat=1 case
		}
		
		srs.el.res=srs.el(srsx,srsy,nx,ny)
		
		## confidence interval
		srs.el.ci=c(srs.el.res$Low,srs.el.res$Up)
		## converage probability
		srs.el.cp=(srs.el.ci[2]>=AUC)*(srs.el.ci[1]<=AUC)
		## length
		srs.el.len=srs.el.res$Up-srs.el.res$Low
		
		# save results
		## confidence interval
		rsse.rnew.ci[sind,]=rss.el.rnew.ci
		rssk.ci[sind,]=rss.ker.ci
		srs.ci[sind,]=srs.el.ci
		## coverage probability
		cp[sind,]=c(rss.el.rnew.cp, rss.ker.cp, srs.el.cp)
		len[sind,]=c(rss.el.rnew.len, rss.ker.len, srs.el.len)
	}
	# return results
	ciavg=rbind(
		apply(rsse.rnew.ci,2,mean),
		apply(rssk.ci,2,mean),
		apply(srs.ci,2,mean)
	)
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}


########################################################
# lognormal
########################################################
main.brss.lognormal=function(pind,params,nsim=10000){
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
		set.seed(30000*pind + sind)
		#############################################################
		# RSS
		rss.cond=TRUE
		while(rss.cond==T){
			# dist of X
			nxsam = 10000
			x=stats::rlnorm(nxsam)
			# concomitant of X
			mux=exp(1/2)
			sdx=sqrt( (exp(1)-1)*exp(1) )
			cx=corx*((x-mux)/sdx) + sqrt(1-corx^2)*rnorm(nxsam)
			
			# dist of Y
			ny.sam= 10000
			sig=1
			muy.temp = sqrt(sig^2+1)*qnorm(AUC)
			sdy.temp = sig
			muy=exp(muy.temp+(sdy.temp^2)/2)
			sdy=sqrt( (exp(sdy.temp^2)-1)*exp(2*muy.temp+sdy.temp^2) )
			y=stats::rlnorm(ny.sam,meanlog=muy.temp,sdlog=sdy.temp)
			# concomitant of y
			cy=cory*((y-muy)/sdy) + sqrt(1-cory^2)*rnorm(ny.sam)
			
			rssx=RSSampling::con.rss(x,cx,m=m,r=k)
			rssy=RSSampling::con.rss(y,cy,m=n,r=l)
			
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
											is.na(rss.el.rnew.res$Low),
											is.na(rss.el.rnew.res$Up),
											is.na(rss.ker.ci[1]),
											is.na(rss.ker.ci[2]))
			)
		}
		
		#############################################################
		# SRS
		srs.cond=TRUE
		while (srs.cond==T) {
			srsx=stats::rlnorm(nx)
			srsy=stats::rlnorm(ny,meanlog=muy.temp,sdlog=sdy.temp)
			
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
		
		# save results
		## confidence interval
		rsse.rnew.ci[sind,]=rss.el.rnew.ci
		rssk.ci[sind,]=rss.ker.ci
		srs.ci[sind,]=srs.el.ci
		## coverage probability
		cp[sind,]=c(rss.el.rnew.cp, rss.ker.cp, srs.el.cp)
		len[sind,]=c(rss.el.rnew.len, rss.ker.len, srs.el.len)
	}
	# return results
	ciavg=rbind(
		apply(rsse.rnew.ci,2,mean),
		apply(rssk.ci,2,mean),
		apply(srs.ci,2,mean)
	)
	cpavg=apply(cp,2,mean)
	lenavg=apply(len,2,mean)
	return(list(ci=ciavg,cp=cpavg,len=lenavg))
}
