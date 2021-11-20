get.mut<-function(dat,tm,k=10,n.pred=NULL){
  gams<-apply(dat,2,function(y,x){gm<-mgcv:::gam(y~s(x,k=k));return(gm)},x=tm)
  if (is.null(n.pred)){
    n.pred=length(tm)^2
  }
  x<-seq(min(tm),max(tm),length=n.pred)
  return(list(mut=lapply(gams,predict,newdata=data.frame(x=x),type="response"),t=x))
}

get.G<-function(mut,tm,VH,p,r){
  G<-NULL
  G[[1]]<-matrix(0,ncol=r,nrow=p)
  VH.o<- VH(mut[1,])
  for (i in 2:length(tm)){
    VH.n<-VH(mut[i,])
    G[[i]]<-G[[i-1]]+(VH.o+VH.n)*(tm[i]-tm[i-1])/2
    VH.o<-VH.n
  }
  return(G)
}

get.GtG<-function(G,tm,r){
  GtG<-matrix(0,ncol=r,nrow=r)
  GtG.o<- t(G[[1]])%*%G[[1]]
  for (i in 2:length(tm)){
    GtG.n<-t(G[[i]])%*%G[[i]]
    GtG<-GtG+(GtG.o+GtG.n)*(tm[i]-tm[i-1])/2
    GtG.o<-GtG.n
  }
  return(GtG)
}

get.Gt<-function(G,tm,p,r){
  Gt<-matrix(0,ncol=p,nrow=r)
  Gt.o<- t(G[[1]])
  for (i in 2:length(tm)){
    Gt.n<-t(G[[i]])
    Gt<-Gt+(Gt.o+Gt.n)*(tm[i]-tm[i-1])/2
    Gt.o<-Gt.n
  }
  return(Gt)
}

get.mutG<-function(mut,G,tm,r){
  mutG<-matrix(0,ncol=1,nrow=r)
  mutG.o<- t(G[[1]])%*%mut[1,]
  for (i in 2:length(tm)){
    mutG.n<-t(G[[i]])%*%mut[i,]
    mutG<-mutG+(mutG.o+mutG.n)*(tm[i]-tm[i-1])/2
    mutG.o<-mutG.n
  }
  return(mutG)
}


get.intmu<-function(mut,tm,p){
  intmu<-rep(0,ncol(mut))
  intmu.o<- mut[1,]
  for (i in 2:length(tm)){
    intmu.n<-mut[i,]
    intmu<-intmu+(intmu.o+intmu.n)*(tm[i]-tm[i-1])/2
    intmu.o<-intmu.n
  }
  return(intmu)
}



###########


gXn = function(X1,X2){
  gXn = matrix(c(X1,0,
                 -X1*X2,
                 X1*X2,
                 0,
                 -X2),
               nrow=2,ncol=3)
  return(gXn)
}

Xn1_funcion = function(llim,rlim,vals){
  dat = vals$times<rlim & vals$times>llim
  vals = vals$X1[dat]
  return(mean(vals))
}
Xn2_funcion = function(llim,rlim,vals){
  dat = vals$times<rlim & vals$times>llim
  vals = vals$X2[dat]
  return(mean(vals))
}

par_estimation <- function(out,I){

  n = I*I
  i = 1:I
  ct = out$t[length(out$t)]
  delta = ct/n
  Delta = ct/I
  Si = i * Delta
  Sim1 = c(0,Si[-length(Si)])

  ti = (1:n) * delta
  vals = NULL
  x1 = as.numeric(out$x[,1])
  x2 = as.numeric(out$x[,2])
  t1 = as.numeric(out$t)
  for(j in 1:(n-1)){
    vals = rbind(vals,
                 data.frame(times=ti[j],
                            X1=x1[min(which(t1>ti[j]))-1],
                            X2=x2[min(which(t1>ti[j]))-1]))
  }
  X1 = X2 = NULL
  Gn_i = gn = list()
  GnT = matrix(0,nrow=2,ncol=3)
  An = matrix(0,nrow=2,ncol=3)
  Bn = matrix(0,nrow=3,ncol=3)
  s1 = matrix(0,nrow=3,ncol=1)
  x12 = c(0,0)
  for(j in i){
    X1[j] = Xn1_funcion(llim=Sim1[j],rlim=Si[j],vals=vals)
    X2[j] = Xn2_funcion(llim=Sim1[j],rlim=Si[j],vals=vals)

    gn[[j]] = gXn(X1[j],X2[j])
    GnT = GnT+gn[[j]]
    Gn_i[[j]] = GnT
    An =  An + Gn_i[[j]]*Delta
    Bn =  Bn + (t(Gn_i[[j]])%*%Gn_i[[j]])*Delta
    s1 = s1 + (t(Gn_i[[j]])%*%c(X1[j],X2[j]))*Delta
    x12 = x12 + c(X1[j],X2[j])

  }

  om1 = rbind(cbind(diag(rep(ct,2)),An),cbind(t(An),Bn))
  om2 = c(x12*Delta,s1)
  omega = solve(om1,om2)

  return(omega)
}

rate<-function(y,K,theta){
  theta*apply(K,2,function(k,n){prod(choose(n,k))},n=y)
}

odest<-function(dat,tms,V,K=NULL,k=10,n.pred=NULL,plot.it=FALSE, QP=TRUE){
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

  mu <- get.mut(dat,tms,k=k)
  mut<- matrix(unlist(mu$mut),ncol=p)
  tm <- mu$t
  intmu<-get.intmu(mut,tm,p = p)
  G <- get.G(mut,tm,VH,p,r)
  GtG <- get.GtG(G,tm,r = r)
  Gt <- get.Gt(G,tm,p,r)
  mutG<- get.mutG(mut,G,tm,r = r)
  dtI <- diag(rep(max(tm)-min(tm),p))
  D <- rbind(cbind(GtG,Gt), cbind(t(Gt),dtI))
  d <- c(mutG,intmu)
  if (!QP){
    pars<- solve(D,d)
  } else {
    fit <- quadprog:::solve.QP(D, d,diag(rep(1,r+p)),rep(0,r+p))
    pars <- fit$solution
  }
  return(list(theta=pars[1:r],mu0=pars[r+1:p],D=D,d=d))
}

odest2<-function(dat,tms,V,K,k=10,n.pred=NULL,plot.it=FALSE, QP=TRUE){
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

  I = length(tms)

  #intmu = matrix(0,nrow = I, ncol = p)
  #mutG = matrix(0,nrow = I, ncol = r)
  D = matrix(0,nrow = p+r, ncol = p+r)
  d = matrix(0,nrow = p+r, ncol = 1)

  #GtG = Gt = list()

  M=m=list()

  mu <- get.mut(dat = dat,tm = tms,k=k)
  mut<- matrix(unlist(mu$mut),ncol=p)
  tms <- mu$t

  G <- get.G(mut,tms,VH,p,r)


  for(j in 1:I){
    tm = tms[1:I+I*(j-1)]
   # da = dat[1:I+I*(j-1),]
    submu = mut[1:I+I*(j-1),]
    subG = G[1:I+I*(j-1)]
    intmu <- get.intmu(mut = submu,tm = tm,p = p)

    GtG <- get.GtG(subG,tm,r = r)
    Gt <- get.Gt(subG,tm,p,r)
    mutG<- get.mutG(submu,subG,tm,r = r)
    dtI <- diag(rep(max(tm)-min(tm),p))
    M[[j]] <- rbind(cbind(GtG,Gt), cbind(t(Gt),dtI))
    m[[j]] <- c(mutG,intmu)

    D = D + M[[j]]
    d = d + m[[j]]
  }

  if (!QP){
    pars<- solve(D,d)
  } else {
    fit <- quadprog:::solve.QP(D, d,diag(rep(1,r+p)),rep(0,r+p))
    pars <- fit$solution
  }
  A = D/I
  B = matrix(0,nrow=p+r,ncol=p+r)
  for(j in 1:I){
    B = B + (M[[j]]%*%pars-m[[j]])%*%t(M[[j]]%*%pars-m[[j]])*(1/(I-1))
  }
  Var = (solve(A)%*%B%*%solve(A))
  return(list(theta=pars[1:r],mu0=pars[r+1:p],Var=Var,D=D,d=d,A=A,B=B,I=I))
}

