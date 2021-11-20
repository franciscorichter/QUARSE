

plot(stepfun(out$t,out$x[,1]),pch="")
plot(stepfun(out$t,out$x[,2]),pch="")
plot(out$x,type="l")


#####################
##  Estimation
#####################


I = 30
out <- gillespie(N,10000)
par_estimation(out,I)

Omega = NULL
for(h in 1:100){
  print(h)

  out <- NULL
  attempt <- 1
  while( is.null(out) && attempt <= 5 ) {
    attempt <- attempt + 1
    try(
      out <- gillespie(N,10000)
    )
  }
  omega = par_estimation(out,I)
  Omega = rbind(Omega,omega)
}




xInitial <- est.qp$mu0
theta <- est.qp$theta
times <- seq(min(tm),max(tm),length=1000)

out   <- ode(xInitial,times,ozone.ode,theta,method='ode45')

I=information_matrix(dat,out)
