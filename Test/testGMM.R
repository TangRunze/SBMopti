rm(list=ls())
setwd("C:/Dropbox/Github/SBMopti/")

require(mclust)

d = 3

asge = function(Ai,d)
  # this function generates the adjacency spectral graph embedding
{
  S = eigen(Ai)
  Xhat = S$vectors[,1:d] %*% diag(sqrt(S$values[1:d]))
}

#300 mclust

MaxIter = 1000
m = 100

set.seed(123)
require(combinat)
Lmclust = rep(0,MaxIter)
classesi = c(rep(1,m), rep(2,m), rep(3,m))
for (i in 1:MaxIter) {
  print(i)
  
  Ai = as.matrix(read.csv(paste("graph",i,sep=""),sep = " ",header = F))
  Ai = Ai[1:300,1:300]

  #Xhat = as.matrix(read.csv(paste("xhat",i,sep=""),sep = " ",header = F))
  #Xhat = Xhat[,1:3]
  #Xhat = matrix(t(Xhat), 300, 3)
  
  Xhat = asge(Ai, d)
  
  M = Mclust(Xhat,d,modelNames=c("EII"))
  perm = permn(1:3)
  L = 3*m
  for (iter in 1:length(perm)) {
    classtmp = c(rep(perm[[iter]][1],m),rep(perm[[iter]][2],m),rep(perm[[iter]][3],m))
    if (mean(M$classification != classtmp)<L) {
      L = mean(M$classification != classtmp)
    }
  }
  Lmclust[i] = L
}
mean(Lmclust)
#95%CI
c(mean(Lmclust)-qnorm(0.95)*sqrt(var(Lmclust))/sqrt(MaxIter),mean(Lmclust)+qnorm(0.95)*sqrt(var(Lmclust))/sqrt(MaxIter))


#300 mclust

MaxIter = 500
m = 100

set.seed(123)
require(combinat)
Lmclust = rep(0,MaxIter)
classesi = c(rep(1,m), rep(2,m), rep(3,m))
for (i in 1:MaxIter) {
  print(i)
  
  Xhat = as.matrix(read.csv(paste("xhat",i,sep=""),sep = " ",header = F))
  Xhat = Xhat[,1:3]
  Xhat = matrix(t(Xhat), 300, 3)
  
  tau = as.matrix(read.csv(paste("tauhat",i,sep=""),sep = " ",header = F))
  tau = tau[,1:300]
  
  perm = permn(1:3)
  L = 3*m
  for (iter in 1:length(perm)) {
    classtmp = c(rep(perm[[iter]][1],m),rep(perm[[iter]][2],m),rep(perm[[iter]][3],m))
    if (mean(tau != classtmp)<L) {
      L = mean(tau != classtmp)
    }
  }
  Lmclust[i] = L
}
mean(Lmclust)
#95%CI
c(mean(Lmclust)-qnorm(0.95)*sqrt(var(Lmclust))/sqrt(MaxIter),mean(Lmclust)+qnorm(0.95)*sqrt(var(Lmclust))/sqrt(MaxIter))


