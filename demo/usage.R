
library(mgcv)
library(quadprog)
library(mvtnorm)
library(deSolve)

library(meanfieldapprox)

V=matrix(c(2,-1,0,-1,-1,1,1,1,-1,-1,2,-1),ncol=4)
K = -V
K[K<0]=0
o=c(200,100,200)

theta=c(1/2,1/4,12,1)

##
V=matrix(c(1,0,-1,1,0,-1),ncol=3)
K=matrix(c(1,0,1,1,0,1),ncol=3)
o=c(200,500)
theta=c(1,0.005,0.6)
##

set.seed(2)
n.sim<-2000
ozone.gil <- react(theta = theta,V =  V,K = K, y0=o,n.iter=2000)

sel<-round(seq(350,405,length=80))
dat<-ozone.gil$state[sel,]
tm<-ozone.gil$tm[sel]

est.raw<-odest2(dat = dat,tms = tm,V = V,K = K,k=5,QP=TRUE)
#est.qp<-odest(dat = dat,tms = tm,V = V,K = K,k=5,QP=TRUE)
est.qp<-odest(dat = dat,tms = tm,V = V,K = K,k=5,QP=FALSE)

xInitial <- est.raw$mu0
theta <- est.raw$theta
times <- seq(min(tm),max(tm),length=1000)
out   <- ode(xInitial,times,lv.ode,theta,method='ode45')

I=information_matrix(dat,out,tm)
det(I)
sum(diag(I))

b=bootstrap_ode(dat = dat,tm = tm,V = V,K = K,nsim = 100)


