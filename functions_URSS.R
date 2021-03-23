library(emplik)
library(tidyverse)
library(RSSampling)
library(snow)

####################################################################################
####### el.rss.test(), from Owen, Modified by Mai Zhou, Modifed by Chul Moon #######
####################################################################################

el.rss.test <- function(ou, mu, lam, maxit=25, gradtol=1e-7, 
												svdtol = 1e-9, itertrace=FALSE ){
	# input is OU, 1-U_{rs}
	x <- as.matrix(unlist(ou))
	xn <- nrow(x)
	xp <- ncol(x)
	mu <- as.vector(mu)
	if( length(mu) !=xp )
		stop("mu must have same dimension as observation vectors.")
	if( xn <= xp )
		stop("Need more observations than length(mu) in el.test().")
	
	sou = list() # standardized OU
	for (ii in 1:length(ou)){
		sou[[ii]]=(ou[[ii]]-mu)/length(ou[[ii]])
	}
	z = as.matrix(unlist(sou))
	
	#
	#    Scale the problem, by a measure of the size of a 
	# typical observation.  Add a tiny quantity to protect
	# against dividing by zero in scaling.  Since z*lam is
	# dimensionless, lam must be scaled inversely to z.
	#
	TINY <- sqrt( .Machine$double.xmin )
	scale <- mean( abs(z) ) + TINY
	z <- z/scale
	if( !missing(lam) ){
		lam <- as.vector(lam)
		lam <- lam*scale
		if( logelr(z,rep(0,xp),lam)>0 )lam <- rep(0,xp)
	}
	if(  missing(lam)  )
		lam <- rep(0,xp)
	#
	#     Take some precaution against users specifying
	# tolerances too small.
	#
	
	if( svdtol < TINY )svdtol <- TINY
	if( gradtol < TINY)gradtol <- TINY
	
	#
	#    Preset the weights for combining Newton and gradient
	# steps at each of 16 inner iterations, starting with
	# the Newton step and progressing towards shorter vectors
	# in the gradient direction.  Most commonly only the Newton
	# step is actually taken, though occasional step reductions
	# do occur.
	#
	
	nwts <- c( 3^-c(0:3), rep(0,12) )
	gwts <- 2^( -c(0:(length(nwts)-1)))
	gwts <- (gwts^2 - nwts^2)^.5
	gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
	
	#
	#    Iterate, finding the Newton and gradient steps, and
	# choosing a step that reduces the objective if possible.
	#
	
	nits <- 0
	gsize <- gradtol + 1
	while(  nits<maxit && gsize > gradtol  ){
		arg  <- 1 + z %*% lam
		wts1 <- as.vector( llogp(arg, 1/xn) )
		wts2 <- as.vector( -llogpp(arg, 1/xn) )^.5
		grad <- as.matrix( -z*wts1 )
		#############grad <- as.vector( apply( grad, 2, sum ) )
		grad <- as.vector(rowsum(grad, rep(1, nrow(grad)) ) )
		gsize <- mean( abs(grad) )
		hess <- z*wts2
		#                                   -1
		#    The Newton step is -(hess'hess)    grad,
		#  where the matrix hess is a sqrt of the Hessian.
		#  Use svd on hess to get a stable solution.
		#
		
		## may try La.svd() in R (v. > 1.0) for better LAPACK.
		## or use QR decomposition on hess to solve it.
		
		svdh <- svd( hess )
		##  svdh <- La.svd( hess )
		if( min(svdh$d) < max(svdh$d)*svdtol )
			svdh$d <- svdh$d + max(svdh$d)*svdtol
		nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
		## nstep <- t(svdh$vt) %*% (t(svdh$u)/svdh$d)
		nstep <- as.vector( nstep %*% matrix(wts1/wts2,xn,1) )
		
		gstep <- -grad
		if(  sum(nstep^2) < sum(gstep^2) )
			gstep <- gstep*(sum(nstep^2)^.5/sum(gstep^2)^.5)
		ologelr <- -sum( llog(arg,1/xn) )
		ninner <- 0
		for(  i in 1:length(nwts) ){
			nlogelr <- logelr( z,rep(0,xp),lam+nwts[i]*nstep+gwts[i]*gstep )
			if( nlogelr < ologelr ){
				lam <- lam+nwts[i]*nstep+gwts[i]*gstep
				ninner <- i
				break
			}
		}
		nits <- nits+1
		if(  ninner==0  )nits <- maxit
		if( itertrace )
			print( c(lam, nlogelr, gsize, ninner) )
	}
	
	list( "-2LLR" = -2*nlogelr, Pval = 1-pchisq(-2*nlogelr, df=xp),
				lambda = lam/scale, grad=grad*scale,
				hess=t(hess)%*%hess*scale^2, wts=wts1, nits=nits )
}

logelr <- function( x, mu, lam ){ 
	x <- as.matrix(x)
	xn <- nrow(x)
	xp <- ncol(x)
	if(  xn <= xp  )
		stop("Need more observations than variables in logelr.")
	mu <- as.vector(mu)
	if(  length(mu) != xp  )
		stop("Length of mean doesn't match number of variables in logelr.")
	
	z <- t( t(x) -mu )
	arg <- 1 + z %*% lam
	return( - sum( llog(arg,1/xn) ) ) 
}

#
#    The function llog() is equal to the natural
#  logarithm on the interval from eps >0 to infinity.
#  Between -infinity and eps, llog() is a quadratic.
#  llogp() and llogpp() are the first two derivatives
#  of llog().  All three functions are continuous
#  across the "knot" at eps.
#
#    A variation with a second knot at a large value
#  M did not appear to work as well.
#
#    The cutoff point, eps, is usually 1/n, where n
#  is the number of observations.  Unless n is extraordinarily
#  large, dividing by eps is not expected to cause numerical
#  difficulty.
#

llog <- function( z, eps ){
	ans <- z
	avoidNA <- !is.na(z)
	lo <- (z<eps) & avoidNA  ### added 3/2012
	ans[ lo  ] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
	ans[ !lo ] <- log( z[!lo] )
	ans
}

llogp <- function( z, eps ){
	ans <- z
	avoidNA <- !is.na(z)    ###added 3/2012
	lo <- (z<eps) & avoidNA
	ans[ lo  ] <- 2.0/eps - z[lo]/eps^2
	ans[ !lo ] <- 1/z[!lo]
	ans
}

llogpp <- function( z, eps ){
	ans <- z
	avoidNA <- !is.na(z) 
	lo <- (z<eps) & avoidNA    ### added same avoidNA as above
	ans[ lo  ] <- -1.0/eps^2
	ans[ !lo ] <- -1.0/z[!lo]^2
	ans
}

##################################################################
##################### functions for RSS + EL #####################
##################################################################

# Lemma
## v
v10fun = function(yvec,xval){
	mean(yvec>=xval)
}
v01fun = function(xvec,yval){
	mean(yval>=xvec)
}

## V10
v10 = function(xij,rssy) {
	mean(unlist(lapply(rssy,v10fun,xij)))
}
## v01
v01 = function(yrs,rssx) {
	mean(unlist(lapply(rssx,v01fun,yrs)))
}

# adjusted EL
el.auc <- function(mu,x,r.adj) {
	elres = el.rss.test(x,mu)
	# adjust log EL
	elres$"-2LLR" = elres$"-2LLR"*r.adj
	return(elres)
}

# unbalanced RSS
urss.el.rnew = function(rssx.vec,rssy.vec,rssx,rssy,m,ki,n,lr,nx,ny) {
	# ki is a vector of size m
	# lr is a vector of size n
	
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
	
	## One minus U
	### for EL
	OU = list()
	for (rr in 1:length(rssy)) {
		OU[[rr]] = apply(as.matrix(rssy[[rr]]),1,v01,rssx)
	}
	
	## sum(1-Uhat-delhat)^2
	### for EL scaling
	OUDsq.temp = rep(NA,length(rssy))
	for (rr in 1:length(rssy)) {
		OUDsq.temp[rr] = mean((apply(as.matrix(rssy[[rr]]),1,v01,rssx)-delhat)^2)
	}
	OUDsq = mean(OUDsq.temp)
	
	## v10_{[i]}.bar
	v10i.bar=function(rssx,rssy,i){
		return(mean(apply(as.matrix(rssx[[i]]),1,v10,rssy)))
	}
	## v01_{[r]}.bar
	v01r.bar=function(rssx,rssy,r){
		return(mean(apply(as.matrix(rssy[[r]]),1,v01,rssx)))
	}
	## (S^10_{[i]})^2
	s10i.sq=function(rssx,rssy,i,ki){
		v10val=apply(as.matrix(rssx[[i]]),1,v10,rssy)
		return(sum( (v10val-v10i.bar(rssx,rssy,i))^2 )/(ki[i]-1) )
	}
	## (S^01_{[r]})^2
	s01r.sq=function(rssx,rssy,r,lr){
		v01val=apply(as.matrix(rssy[[r]]),1,v01,rssx)
		return(sum( (v01val-v01r.bar(rssx,rssy,r))^2 )/(lr[r]-1) )
	}
	## (S^10)^2
	s10.sq.vec = rep(NA,m)
	for (ii in 1:m) {
		s10.sq.vec[ii] = s10i.sq(rssx,rssy,ii,ki)
	}
	s10.sq = mean(s10.sq.vec)
	## (S^01)^2
	s01.sq.vec = rep(NA,n)
	for (rr in 1:n) {
		s01.sq.vec[rr] = s01r.sq(rssx,rssy,rr,lr)
	}
	s01.sq = mean(s01.sq.vec)
	
	## S^2
	S.sq = (ny*s10.sq+nx*s01.sq)/(nx+ny)
	
	if (S.sq==0) return(x=list(Low=NA,Up=NA)) 
	
	else {
		# Theorem 2: adjusted asymptotic dist
		## adjustment
		r.adj = (nx)/(ny+nx)*OUDsq/(S.sq)
		
		rss.el.res = emplik::findUL2(step=0.005, fun=el.auc, MLE=delhat, x=OU, r.adj=r.adj)
		return(rss.el.res)
	}
}