
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
			# RSS sampling of X
		  rssx=con.rss.sim.normal(m=m,r=k,mu=0,sigma=1,cor=corx)
			
			# RSS sampling of Y
			mu = sqrt(5)*qnorm(AUC)
			sigma = 2
			rssy=con.rss.sim.normal(m=n,r=l,mu=mu,sigma=sigma,cor=cory)

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
		  
		  # RSS sampling of X
		  rssx=con.rss.sim.uniform(m=m,r=k,mu=0.5,sigma=sqrt(1/12),cor=corx,theta=1)
		  
		  # RSS sampling of Y
			theta = 1/(2-2*AUC)
			muy = theta/2
			sdy = sqrt(1/12*theta^2)
			
			rssy=con.rss.sim.uniform(m=n,r=l,mu=muy,sigma=sqrt(1/12*theta^2),
			                         cor=cory,theta=theta)

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
			# RSS sampling of X
			mux=exp(1/2)
			sdx=sqrt( (exp(1)-1)*exp(1) )
			
			# RSS sampling of Y
			sig=1
			muy.temp = sqrt(sig^2+1)*qnorm(AUC)
			sdy.temp = sig
			muy=exp(muy.temp+(sdy.temp^2)/2)
			sdy=sqrt( (exp(sdy.temp^2)-1)*exp(2*muy.temp+sdy.temp^2) )

			rssx=con.rss.sim.lognormal(m=m,r=k,meanlog=0,sdlog=1,
			                           mu=mux,sigma=sdx,cor=corx)
			rssy=con.rss.sim.lognormal(m=n,r=l,meanlog=muy.temp,sdlog=sdy.temp,
			                           mu=muy,sigma=sdy,cor=cory)
			
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


########################################################
# function for different reference functions
########################################################
main.brss.normal.reference=function(pind,params,nsim=10000){
  source("./functions_BRSS.R")
  param=params[pind,]
  m=param$m
  n=param$n
  k=param$nx/param$m
  l=param$ny/param$n
  nx=param$nx
  ny=param$ny
  AUC=param$AUC
  corx=cory=param$corxy
  
  len = cp = matrix(NA,nsim,2)
  rsse.rnew.ci = rsse.switch.ci = matrix(NA,nsim,2)
  
  for (sind in 1:nsim) {
    set.seed(23500*pind + sind)
    print(sind)
    #############################################################
    # RSS
    rss.cond=TRUE
    while(rss.cond==T){
      
      rssx = con.rss.sim.normal(m=m,r=k,mu=0,sigma=1,cor=corx)
      rssy = con.rss.sim.normal(m=n,r=l,mu=sqrt(5)*qnorm(AUC),sigma=2,cor=cory)
      
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
      # RSS + EL: F as the reference
      ######################################################
      rss.el.rnew.res=rss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
      ## compute confidence interval
      rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
      ## compute coverage probability
      rss.el.rnew.cp=(rss.el.rnew.ci[2]>=AUC)*(rss.el.rnew.ci[1]<=AUC)
      
      ## save confidence interval
      rsse.rnew.ci[sind,]=rss.el.rnew.ci
      
      ######################################################
      # RSS + EL: G as the reference
      ######################################################
      rss.el.switch.res=rss.el.switch(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
      ## confidence interval
      rss.el.switch.ci=c(rss.el.switch.res$Low,rss.el.switch.res$Up)
      ## coverage probability
      rss.el.switch.cp=(rss.el.switch.ci[2]>=AUC)*(rss.el.switch.ci[1]<=AUC)
      
      ## save confidence interval
      rsse.switch.ci[sind,]=rss.el.switch.ci
      
      ## save coverage probability
      cp[sind,]=c(rss.el.rnew.cp, rss.el.switch.cp)
      
      rss.cond=any( c(delhat==1, # remove delta-hat=1 case
                      is.na(rss.el.rnew.res$Low),
                      is.na(rss.el.rnew.res$Up),
                      is.na(rss.el.switch.res$Low),
                      is.na(rss.el.switch.res$Up))
      )
    }
  }
  # return results
  ciavg=rbind(
    apply(rsse.rnew.ci,2,mean),
    apply(rsse.switch.ci,2,mean)
  )
  
  cilensd=rbind(
    sd(rsse.rnew.ci[,2]-rsse.rnew.ci[,1]),
    sd(rsse.switch.ci[,2]-rsse.switch.ci[,1])
  )
  
  cilenmean=rbind(
    mean(rsse.rnew.ci[,2]-rsse.rnew.ci[,1]),
    mean(rsse.switch.ci[,2]-rsse.switch.ci[,1])
  )
  cpavg=apply(cp,2,mean)
  
  return(list(ci=ciavg,cp=cpavg,lensd=cilensd,lenmean=cilenmean))
}
