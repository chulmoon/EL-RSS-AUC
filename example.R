###############################################
# BRSS example
###############################################
source("./functions_BRSS.R")
# parameters
nx=ny=40 # sample size
m=n=2 # set size
k=l=nx/m # number of cycles
corx=cory=1 # correlation
AUC=0.9 # AUC

set.seed(9)

## dist of X
nxsam = 10000
x=stats::rnorm(nxsam)
## concomitant of X
cx=corx*x + sqrt(1-corx^2)*rnorm(nxsam)
## dist of Y
ny.sam= 10000
mu = sqrt(5)*qnorm(AUC)
sigma = 2
y=stats::rnorm(ny.sam,mu,sigma)
# concomitant of y
cy=cory*((y-mu)/sigma) + sqrt(1-cory^2)*rnorm(ny.sam)

## BRSS sampling
rssx=RSSampling::con.rss(x,cx,m=m,r=k)
rssy=RSSampling::con.rss(y,cy,m=n,r=l)

## balanced ranked set samples without ranking information
rssx.vec=as.vector(rssx$sample.x)
rssy.vec=as.vector(rssy$sample.x)

## delta-hat
delhat = 0
for (i in 1:length(rssx.vec)){
  for (j in 1:length(rssy.vec)) {
    delhat=delhat+(rssy.vec[j]>=rssx.vec[i])
  }
}
delhat = delhat/(m*k*n*l)

## BRSS-EL
rss.el.rnew.res=rss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,k,n,l)
## confidence interval
rss.el.rnew.ci=c(rss.el.rnew.res$Low,rss.el.rnew.res$Up)
print(paste0("Confidence interval of BRSS-EL: (",
             round(rss.el.rnew.ci[1],3),",",
             round(rss.el.rnew.ci[2],3),")"))

###############################################
# URSS example
###############################################
source("./functions_URSS.R")
# parameters
nx=ny=40 # sample size
m=n=2 # set size
prop=0.6 # 60% samples for the second rank
ki=c(nx*0.5,nx*0.5) # number of cycles of X of URSS
lr=c(ny*prop,ny*(1-prop)) # number of cycles of Y of URSS
corx=cory=1 # correlation
AUC=0.9 # AUC

set.seed(9)

## dist of X
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

# URSS sampling
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

# unbalanced ranked set samples without ranking information
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
  
## URSS-EL
urss.el.rnew.res=urss.el.rnew(rssx.vec,rssy.vec,rssx,rssy,m,ki,n,lr,nx,ny)
## confidence interval
urss.el.rnew.ci=c(urss.el.rnew.res$Low,urss.el.rnew.res$Up)
print(paste0("Confidence interval of URSS-EL: (",
             round(urss.el.rnew.ci[1],3),",",
             round(urss.el.rnew.ci[2],3),")"))