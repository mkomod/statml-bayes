#replicate_1_parrallel

#Reproducing results from [1]
setwd('~/phd/project1/')
library('foreach')
library('doParallel')
source('replicate_1_functions.R')
numCores = 20

cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl)
foreach(c = 1:numCores) %dopar% {
  output_path = 'output/'
  source('replicate_1_functions.R')
  library('rmutil')
  library('mvtnorm')
  
  N = c(100,178,316,562,1000)
  K = 3
  
  for(k in 1:K){
    #print(paste0('n = ',n))
    A = matrix(runif(16, 0,10), nrow = 4)
    covNorm = A%*%t(A)
    meanNorm = rep(0,4)
    Bstar = solve(A)  #inverse of A
    cstar = log(abs(det(Bstar))) - log(4)
    thetastar = c(as.vector(matrix(Bstar,ncol=1)),cstar)
    
    thetaflatM = matrix(0,nrow=5,ncol=17 )
    for(i in 1:5){
      A_init = matrix(runif(16, 0,10), nrow = 4)
      B_init = solve(A)
      c_init = log(abs(det(B_init))) - log(4)
      theta_init = c(as.vector(matrix(B_init,ncol=1)),c_init)
      thetaflatM[i,] = theta_init
    }
    
    for(n in N){
      sample_size = n
      local_opt_B = NULL
      local_opt_c = NULL
      s = matrix(rlaplace(sample_size*4,m = 0,s = sqrt(1/2)),nrow=4)
      x = A%*%s
      y= t(rmvnorm(sample_size, meanNorm,covNorm))
      
      for(j in 1:nrow(thetaFlatM)){
        thetaFlat = thetaFlatM[j,]
        optimised = optim(par=thetaFlat, x=x,y=y,Tee=n,fn=objective_function2,method='BFGS')
        MSE_B = mean((optimised$par[1:16]-thetastar[1:16])^2)
        MSE_c = (optimised$par[17]-thetastar[17])^2
        
        local_opt_B = min(local_opt_B, MSE_B)
        local_opt_c = min(local_opt_c, MSE_c)
      }
      write.table(local_opt_B, paste0(output_path,'B_',n,'.csv'), append=TRUE, row.names = F, col.names=F, sep=',')
      write.table(local_opt_c, paste0(output_path,'c_',n,'.csv'), append=TRUE, row.names = F, col.names=F, sep=',')
    }
  }
}
parallel::stopCluster(cl)

