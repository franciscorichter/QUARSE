
a_i=
b_i=
c_i=
d_i=


### Information matrix for MLE variance
ozone.ode <- function(times,x,theta)
{
  dx = x;
  dimnames(dx) = dimnames(x);
  dx[1] = 2*theta[1]*x[2]-theta[2]*x[1]*x[2] + theta[3]*x[3] - theta[4]*x[1]*x[3]
  dx[2] = -theta[1]*x[2] -theta[2]*x[1]*x[2] + theta[3]*x[3] + 2*theta[4]*x[1]*x[3]
  dx[3] = theta[2]*x[1]*x[2] - theta[3]*x[3] - theta[4]*x[1]*x[3]

  return(list(dx))
}

lv.ode <- function(times,x,theta)
{
  dx = x;
  dimnames(dx) = dimnames(x);
  dx[1] = theta[1]*x[1]-theta[2]*x[1]*x[2]
  dx[2] = theta[2]*x[1]*x[2]  - theta[3]*x[2]

  return(list(dx))
}

TF.ode <- function(times,x,theta)
{
  dx = x;
  dimnames(dx) = dimnames(x);
  dx[1] = -theta[1]*x[1]*x[3]                      - theta[3]*x[1]*x[5]                      + theta[5]*x[4]                 + theta[7]*x[6]
  dx[2] =                     - theta[2]*x[2]*x[3]                      - theta[4]*x[2]*x[4]                 + theta[6]*x[5]                + theta[8]*x[6]
  dx[3] = -theta[1]*x[1]*x[3] - theta[2]*x[2]*x[3]                                           + theta[5]*x[4] + theta[6]*x[5]
  dx[4] =  theta[1]*x[1]*x[3]                                           - theta[4]*x[2]*x[4] - theta[5]*x[4]                                + theta[8]*x[6]
  dx[5] =                       theta[2]*x[1]*x[2] - theta[3]*x[1]*x[5]                                      - theta[6]*x[5] + theta[7]*x[6]
  dx[6] =                                            theta[3]*x[1]*x[5] + theta[4]*x[2]*x[4]                                 - theta[7]*x[6] - theta[8]*x[6]

  return(list(dx))
}

information_matrix <- function(dat,out,tm){

  Sigma<-diag(apply(rbind(dat,out[,-1]),2,function(x,n){sum((x[1:n]-x[-(1:n)])^2)/(n-1)},n=nrow(dat)))

  l=1000
  delta = (max(tm)-min(tm))/l

  out=as.data.frame(out)
  names(out) = c("time","mu1","mu2","mu3")

  J1=J2=J3=J4=list()
  for(i in 1:l){
    J1[[i]] = sum(out$mu2[1:i])*delta
    J2[[i]] = sum(out$mu1[1:i]*out$mu2[1:i])*delta
    J3[[i]] = sum(out$mu3[1:i])*delta
    J4[[i]] = sum(out$mu1[1:i]*out$mu3[1:i])*delta
  }

  I = 0
  for(i in 1:l){
    I = I + t(diag(c(J1[[i]],J2[[i]],J3[[i]],J4[[i]])))%*%t(V)%*%solve(Sigma)%*%V%*%diag(c(J1[[i]],J2[[i]],J3[[i]],J4[[i]]))
  }

  return(I)
}
