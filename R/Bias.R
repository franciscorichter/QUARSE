bootstrap_ode <- function(dat,tm,V,K,nsim = 100){

  R = list()
  P = In = NULL
  for(i in 1:nsim){
    #print(i)
    set.seed(i)
    samp = sample(1:nrow(dat),size = nrow(dat),replace = T)
    dat2 = dat[samp,]
    tm2 = tm[samp]
    estimations = odest2(dat = dat2,tms = tm2,V = V,K = K,k=10,QP=T)
    P = rbind(P,estimations$theta)
    In = rbind(In,estimations$mu0)

  }
  bias_theta = colMeans(P)
  bias_initial = colMeans(In)
  var_theta = matrixStats:::colVars(P)
  var_initial = matrixStats:::colVars(In)
  return(list(bias_theta = bias_theta, bias_initial = bias_initial,var_theta=var_theta,var_initial=var_initial,Pars=P,Initial=In))
}


