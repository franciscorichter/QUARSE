rm(list = ls())
# model
V=matrix(c(1,0,-1,1,0,-1),ncol=3)
K=matrix(c(1,0,1,1,0,1),ncol=3)

# parameters
o=c(200,500)
theta=c(1,0.005,0.6)

# syntetic data
#set.seed(3)
n.sim<-2000
ozone.gil <-react(theta = theta,V =  V,K = K, y0=o,n.iter=20000)
#sel<-round(seq(350,405,length=80))
sel=1:10000
dat<-ozone.gil$state[sel,]
tm<-ozone.gil$tm[sel]

plot(dat)

##
N=list()
N$M= c (x1=200,x2=500)
N$Pre= matrix( c(1,0,1,1,0,1), ncol =2,byrow=TRUE)
N$Post= matrix(c(2,0,0,2,0,0), ncol =2,byrow=TRUE)
N$h = function(x,t,th= c(th1=1,th2=0.005,th3=0.6)) {
  with(as.list(c(x,th)),{
    return(c(th1*x1, th2*x1*x2, th3*x2 ))
  })
}
# simulate a realisation of the process and plot it
out = gillespie(N,10000)

dat = out$x
tm = c(0,out$t)

test_AB(V,K,o,theta,dat,tm)

test_AB <- function(V,K,o,theta,dat,tm){

  # number of states
  p = nrow(V)
  # number of reactions
  r = ncol(V)
  VH<-function(y){
    H<-NULL
    for (i in 1:r){
      H<-c(H,prod(choose(y,K[,i])))
    }
    return(V%*%diag(H))
  }
  mu <- get.mut(dat,tm,k=5,n.pred=NULL,plot.it=FALSE)
  mut <- matrix(unlist(mu$mut),ncol=p)
  tmo = tm
  tm <- mu$t

#  G <- get.G(mut,tm,VH,p,r)
  G<-NULL
  G[[1]]<-matrix(0,ncol=r,nrow=p)
  VH.o<- VH(mut[1,])
  #B
  GtG<-matrix(0,ncol=r,nrow=r)
  GtG.o<- t(G[[1]])%*%G[[1]]
  #A
  Gt<-matrix(0,ncol=p,nrow=r)
  Gt.o<- t(G[[1]])

  GnT = matrix(0,nrow=2,ncol=3)
  An = matrix(0,nrow=2,ncol=3)
  Bn = matrix(0,nrow=3,ncol=3)
  s1 = matrix(0,nrow=3,ncol=1)
  x12 = c(0,0)
  for (i in 2:length(tm)){
    VH.n<-VH(mut[i,])
    G[[i]]<-G[[i-1]]+VH.n*(tm[i]-tm[i-1])
    VH.o<-VH.n
    #B
    GtG.n<-t(G[[i]])%*%G[[i]]
    GtG<-GtG+GtG.n*(tm[i]-tm[i-1])

    #A
    Gt.n<-t(G[[i]])
    Gt<-Gt+Gt.n*(tm[i]-tm[i-1])

    GnT = GnT+gXn(mut[i,1],mut[i,2])*(tm[i]-tm[i-1])
    An =  An + GnT*(tm[i]-tm[i-1])
    Bn =  Bn + (t(GnT)%*%GnT)*(tm[i]-tm[i-1])
    s1 = s1 + (t(GnT)%*%c(mut[i,1],mut[i,2]))*(tm[i]-tm[i-1])
    x12 = x12 + c(mut[i,1],mut[i,2])*(tm[i]-tm[i-1])
  }

  B = GtG
  A = Gt

  t2 = all.equal(B,Bn)
  t1 = all.equal(t(A),An)

  dtI <- diag(rep(max(tm)-min(tm),p))

  om1 = rbind(cbind(dtI,An),cbind(t(An),Bn))
  D <- rbind(cbind(B,A), cbind(t(A),dtI))

  om2 = c(x12,s1)
  omega = solve(om1,om2)

  mutG<- get.mutG(mut,G,tm,p,r)
  intmu<-get.intmu(mut,tm,p,r)
  d <- c(mutG,intmu)

  fit <- quadprog:::solve.QP(D, d,diag(rep(1,r+p)),rep(0,r+p))
  pars <- fit$solution

  pars<- solve(D,d)

  return(list(test1=t1,test2=t2,A=A,B=B))
}


