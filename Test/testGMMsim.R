rm(list=ls())
setwd("E:/Dropbox/Github/SBMopti/")
#setwd("C:/Dropbox/Github/SBMopti/")
#setwd("~/Dropbox/GitHub/SBMopti")

require(mclust)

MaxIter = 500
n = 150
eps = 0.1
d = 3

asge = function(Ai,d)
  # this function generates the adjacency spectral graph embedding
{
  S = eigen(Ai)
  Xhat = S$vectors[,1:d] %*% diag(sqrt(S$values[1:d]))
}

# R Mclust
m = n/d

set.seed(123)
require(combinat)
Lmclust = rep(0,MaxIter)
classesi = c(rep(1,m), rep(2,m), rep(3,m))

#taudiff = rep(0, MaxIter)

for (i in 1:MaxIter) {
  print(i)
  
  #Ai = as.matrix(read.csv(paste("./testdata/sim-graph",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
  #Ai <- matrix(Ai[1:n,1:n], n, n)
  
  Xhat = as.matrix(read.csv(paste("./testdata/sim-xhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
  Xhat = array(Xhat, c(n, d))
  
  #Xhat = asge(Ai, d)
  
  #M = Mclust(Xhat,d,modelNames=c("VVV"))
  #M = Mclust(Xhat,d,modelNames=c("VVV"),control=emControl(tol = c(1-1e-16,1-1e-16)))
  M = Mclust(Xhat,d,modelNames=c("VVV"),control=emControl(itmax = 1))
  #hcTree <- hc(modelName = "VVV", data = Xhat)
  #tauinit <- hclass(hcTree, 3)
  #write(tauinit, paste("./testdata/sim-tauinit",i,"-n",n,"-eps",eps,sep=""), ncol=1)
  write(M$classification, paste("./testdata/sim-tauinit",i,"-n",n,"-eps",eps,sep=""), ncol=1)
  
  #taudiff[i] = sum(tauinit != M$classification)
  
  M = Mclust(Xhat,d,modelNames=c("VVV"))
  
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


# Matlab fitgmdist
m = n/d
require(combinat)
Lmclust = rep(0,MaxIter)
classesi = c(rep(1,m), rep(2,m), rep(3,m))
for (i in 1:MaxIter) {
  print(i)
  
  Xhat = as.matrix(read.csv(paste("./testdata/sim-xhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
  Xhat = array(Xhat, c(n, d))
  
  tau = as.matrix(read.csv(paste("./testdata/sim-tauhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
  tau = array(tau, n)
  
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
