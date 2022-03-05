library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)

##################################################################
#################### functions for BRSS + EL #####################
##################################################################

# Lemma
## V10
v10 = function(xij,y,n,l) {
	return(1/(n*l)*sum(xij<=y))
}
v01 = function(yrs,x,m,k) {
	return(1/(m*k)*sum(yrs>=x))
}

# adjusted EL
el.auc <- function(mu,x,r.adj) {
	elres = emplik::el.test(x,mu)
	# adjust log EL
	elres$"-2LLR" = elres$"-2LLR"*r.adj
	return(elres)
}

rss.el.rnew = function(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l) {
  # rss.el.rnew finds the Wilks confidence interval using emplik::findUL2.
  # The adjusted EL function el.auc is used. 
  
  ## input:
  ### rssx and rssy are the ranked set samples obtained by RSSampling::con.rss of x and y, respectively.
  ### rssx.vec and rssy.vec are the sampled observations of rssx and rssy without rank information.
  ### m and n are the set sizes for x and y, respectively.
  ### l and r are the number of cycles for x and y, respectively.
  ## output:
  ### confidence interval obtained by emplik::findUL2
  
	## delta hat
	delhat = 0
	for (i in 1:length(rssx.vec)){
		for (j in 1:length(rssy.vec)) {
			delhat=delhat+(rssy.vec[j]>=rssx.vec[i])
		}
	}
	delhat = delhat/(m*k*n*l)
	
	## One minus U
	OU = rep(NA,length(rssy.vec))
	for (j in 1:length(rssy.vec)) {
		OU[j] = (sum(rssy.vec[j]>=rssx.vec)/(length(rssx.vec)))
	}
	
	## sum(1-Uhat-delhat)^2
	OUDsq = 0
	for (j in 1:length(rssy.vec)) {
		OUDsq = OUDsq + (sum(rssy.vec[j]>=rssx.vec)/(length(rssx.vec))-delhat)^2
	}
	
	v10i.bar=function(rssx,rssy,i){
		return(mean(apply(as.matrix(rssx$sample.x[,i]),1,v10,y=rssy.vec,n=n,l=l)))
	}
	v01r.bar=function(rssx,rssy,r){
		return(mean(apply(as.matrix(rssy$sample.x[,r]),1,v01,x=rssx.vec,m=m,k=k)))
	}
	
	## S^2
	s10i.sq=function(rssx,rssy,i,k){
		v10val=apply(as.matrix(rssx$sample.x[,i]),1,v10,y=rssy.vec,n=n,l=l)
		return(sum( (v10val-v10i.bar(rssx,rssy,i))^2 )/(k-1) )
	}
	s01r.sq=function(rssx,rssy,r,l){
		v01val=apply(as.matrix(rssy$sample.x[,r]),1,v01,x=rssx.vec,m=m,k=k)
		return(sum( (v01val-v01r.bar(rssx,rssy,r))^2 )/(l-1) )
	}
	
	s10.sq.vec = rep(NA,m)
	for (ii in 1:m) {
		s10.sq.vec[ii] = s10i.sq(rssx,rssy,ii,k)
	}
	s10.sq = mean(s10.sq.vec)
	
	s01.sq.vec = rep(NA,n)
	for (rr in 1:n) {
		s01.sq.vec[rr] = s01r.sq(rssx,rssy,rr,l)
	}
	s01.sq = mean(s01.sq.vec)
	
	S.sq = (n*l*s10.sq+m*k*s01.sq)/(m*k+n*l)
	
	if (S.sq==0) return(x=list(Low=NA,Up=NA)) 
	
	else {
		# Theorem 1: adjusted asymptotic dist
		## adjustment
		r.adj = (m*k)/(m*k+n*l)*OUDsq/(n*l*S.sq)
		
		rss.el.res = emplik::findUL2(step=0.005, fun=el.auc, MLE=delhat, x=OU, r.adj=r.adj)
		return(rss.el.res)
	}
}



##################################################################
################## functions for BRSS + KERNEL ###################
##################################################################

## bandwidth
bd = function(x) {
	h=0.9*min(sd(x),IQR(x)/1.34)*length(x)^(-0.2)
	return(h)
}

# AUC-hat
rssker = function(x,y){
	hx = bd(x)
	hy = bd(y)
	auchat = 0
	for (i in 1:length(x)) {
		for (j in 1:length(y)) {
			auchat = auchat + pnorm( (y[j]-x[i])/sqrt(hx^2+hy^2) )
		}
	}
	auchat = auchat/length(x)/length(y)
	return(auchat)
}

### Dr(X_i)
drx = function(xij,yr,hx,hy) {
	dr = 0
	for (j in 1:length(yr)){
		dr = dr + pnorm( (yr[j]-xij)/sqrt(hx^2+hy^2) )
	}
	dr = dr/length(yr)
	return(dr)
}

### Di(Y_r)
diy = function(yrs,xi,hx,hy) {
	di = 0
	for (i in 1:length(xi)){
		di = di + pnorm( (yrs-xi[i])/sqrt(hx^2+hy^2) )
	}
	di = di/length(xi)
	return(di)
}

logit = function(p) log(p/(1-p))

rss.ker = function(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l){
	hx = bd(rssx.vec)
	hy = bd(rssy.vec)
	
	aucker = rssker(rssx.vec,rssy.vec)
	
	# sigma_tilda
	## sigma.110
	### sum of Dr(X_i)
	sumdrx = rep(NA,n)
	for (nn in 1:n){
		sumdrx[nn] = 
			sum( ( purrr::map_dbl(rssx.vec,drx,rssy$sample.x[,nn],hx,hy) - # Ds(X_i)'s
						 	sum( purrr::map_dbl(rssx.vec,drx,rssy$sample.x[,nn],hx,hy))/length(rssx.vec) )^2 )/
			(length(rssx.vec)-1)
	}
	sigma.110 = sum(sumdrx)/(n^2)
	## sigma.011
	### Sum of Di(Y_r)
	sumdiy = rep(NA,m)
	for (mm in 1:m){
		sumdiy[mm] = 
			sum( ( purrr::map_dbl(rssy.vec,diy,rssx$sample.x[,mm],hx,hy) - # Ds(X_i)'s
						 	sum( purrr::map_dbl(rssy.vec,diy,rssx$sample.x[,mm],hx,hy))/length(rssy.vec) )^2 )/
			(length(rssy.vec)-1)
	}
	sigma.011 = sum(sumdiy)/(m^2)
	
	## sigma.tilda
	sigma.tilda = sigma.110/(m*k) + sigma.011/(n*l)
	
	### CI computation
	### simple calculation
	simple.UL = min(aucker + qnorm(0.975)*sqrt(sigma.tilda),1)
	simple.LL = max(aucker - qnorm(0.975)*sqrt(sigma.tilda),0)
	rss.ker.res = c( simple.LL, simple.UL )
	### logit transformation
	#LL=logit(aucker)-qnorm(0.975)*sqrt(sigma.tilda)/(aucker)/(1-aucker)
	#UL=logit(aucker)+qnorm(0.975)*sqrt(sigma.tilda)/(aucker)/(1-aucker)
	#rss.ker.res = c( exp(LL)/(1+exp(LL)), exp(UL)/(1+exp(UL)) )
	return(rss.ker.res)
}

##################################################################
#################### functions for SRS + EL ######################
##################################################################

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


