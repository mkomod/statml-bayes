#Reproducing results from [1]
setwd('~/phd/project1/')
library('rmutil')
library('mvtnorm')
source('replicate_1_functions.R')

#The sampled data
set.seed(123)
sample_size = 100
s = matrix(rlaplace(sample_size*4,m = 0,s = sqrt(1/2)),nrow=4)
A = matrix(runif(16, -10,10), nrow = 4)
x = A%*%s

#These are the answers we are trying to approximate
Bstar = solve(A)  #inverse of A
cstar = log(abs(det(Bstar)))-log(4)
thetastar = c(as.vector(matrix(Bstar,ncol=1)),cstar)

  
#now generate the random data
covNorm = A%*%t(A)
meanNorm = rep(0,4)
y= t(rmvnorm(sample_size, meanNorm,covNorm))

#init params
thetaFlat = c(runif(16,-1,1),-9)

optimised = optim(par=thetaFlat, x=x,y=y,Tee=sample_size,fn=objective_function2,method='CG')

