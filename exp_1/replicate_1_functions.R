#functions
lnpm <- function(b,c,u){
  return( -sqrt(2)*(  sum(abs(b%*%u) )) + c)
}
lnpn = function(u){
  return(log(dmvnorm(u,meanNorm, covNorm)))
}
G = function(u, theta){
  return(lnpm(theta$B,theta$c,u) - lnpn(u))
}
h = function(u,theta){
  return(1/(1+exp(-G(u,theta))))
}
objective_function = function(thetaFlat,Tee,x,y){
  theta = list('B'=matrix(thetaFlat[1:16], nrow = 4),'c'=thetaFlat[17])
  out = 0
  for(t in 1:Tee){
    out = out + log(h(x[,t], theta)) + log(1-h(y[,t],theta))
  }
  return(-(1/2*Tee)*out )
}

#Version 2
objective_function2 = function(thetaFlat,Tee,x,y){
  out = 0
  for(t in 1:Tee){
    lnpmx = -sqrt(2)*(sum(abs(matrix(thetaFlat[1:16],nrow=4)%*%x[,t]))) + thetaFlat[17]
    lnpmy = -sqrt(2)*(sum(abs(matrix(thetaFlat[1:16],nrow=4)%*%y[,t]))) + thetaFlat[17]
    lnpnx = log(dmvnorm(x[,t],meanNorm,covNorm))
    lnpny = log(dmvnorm(y[,t],meanNorm,covNorm))
    Gy = lnpmy - lnpny
    Gx = lnpmx - lnpnx
    hy = 1/(1+exp(-Gy))
    hx = 1/(1+exp(-Gx))
    out = out + log(hx) + log(1-hy)
  }
  return(-(1/2*Tee)*out )
}

