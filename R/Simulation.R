gillespie <- function(N, n, ...) {
  tt = 0
  x = N$M
  S = t(N$Post-N$Pre)
  u = nrow(S)
  v = ncol(S)
  tvec = vector("numeric",n)
  xmat = matrix(ncol=u,nrow=n+1)
  xmat[1,] = x
  for (i in 1:n) {
    h = N$h(x,tt, ...)
    tt = tt+rexp(1,sum(h))
    j = sample(v,1,prob=h)
    x = x+S[,j]
    tvec[i] = tt
    xmat[i+1,] = x
  }
  return(list(t=tvec, x=xmat))
}


# react<-function(theta,V,K,y0=NULL,n.iter=1000){
#   p = nrow(V)
#   r = ncol(V)
#   tm = 0
#   if (is.null(y0)){
#     state = matrix(rpois(p,theta+1),nrow=1)
#   } else {
#     state = matrix(y0,nrow=1)
#   }
#   i = 2
#   while (i <= n.iter){
#     h <- rate(state[i-1,],K,theta)
#     if (sum(h)>0){
#       tm.inc <- rexp(1,sum(h))
#       tm <- c(tm,tm[i-1]+tm.inc)
#       j <- which(rmultinom(1,1,h)==1)
#       state <- rbind(state,state[i-1,]+V[,j])
#       i = i+1
#     } else {
#       i = n.iter+1
#     }
#   }
#   rng<-range(state)
#   return(list(state=state,tm=tm))
# }

react<-function(theta,V,K=NULL,y0=NULL,method=c("ssa","approx"),n.iter=1000,dt=1,q=NULL){
  if (is.null(K)){
    K = -V
    K[K<0]=0
  }
  p = nrow(V)
  r = ncol(V)
  tm = 0
  if (is.null(y0)){
    state = matrix(rpois(p,theta+1),nrow=1)
  } else {
    state = matrix(y0,nrow=1)
  }
  i = 2
  if (method[1]=="ssa"){
    #print("Perform Stochastic Simulation Algorithm")
    while (i <= n.iter){
      h <- c(rate(state[i-1,],K,theta),q)
      if (sum(h)>0){
        tm.inc <- rexp(1,sum(h))
        tm <- c(tm,tm[i-1]+tm.inc)
        j <- which(rmultinom(1,1,h)==1)
        state <- rbind(state,state[i-1,]+V[,j])
        i = i+1
      } else {
        i = n.iter+1
      }
    }
  } else {
    # Euler-Maruyama approximation
   # print("Performing Euler-Maruyama approximation")
    while (i <= n.iter){
      h <-  c(rate(state[i-1,],K,theta),q)
      if (sum(h)>0){
        tm.inc <- dt
        tm <- c(tm,tm[i-1]+tm.inc)
        mu <- tm.inc*V%*%matrix(h,ncol=1)
        beta<- tm.inc*V%*%diag(h)%*%t(V)
        newstate <- state[i-1,]+rmvnorm(1,mu,beta)
        newstate <- pmax(newstate,0)
        state <- rbind(state,newstate)
        i = i+1
      } else {
        i = n.iter+1
      }
    }
  }
  #rng<-range(state)
 # plot(tm,state[,1],type="l",ylab="States",xlab="Time",ylim=rng)
  #apply(rbind(2:p,state[,-1]),2,function(y,tm){lines(tm,y[-1],col=y[1])},tm=tm)
  return(list(state=state,tm=tm))
}

