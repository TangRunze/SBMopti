rm(list=ls())
#setwd("E:/Dropbox/Github/SBMopti/")
setwd("C:/Dropbox/Github/SBMopti/")

require(mclust)

i = 1
n = 150
eps = 0.1
d = 3

m = n/d

set.seed(123)
require(combinat)
classesi = c(rep(1,m), rep(2,m), rep(3,m))

Xhat = as.matrix(read.csv(paste("./testdata/sim-xhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
Xhat = array(Xhat, c(n, d))

M = Mclust(Xhat,d,modelNames=c("VVV"))

hcTree <- hc(modelName = "VVV", data = Xhat)
tauinit <- hclass(hcTree, 3)

perm = permn(1:3)
L = 3*m
for (iter in 1:length(perm)) {
  classtmp = c(rep(perm[[iter]][1],m),rep(perm[[iter]][2],m),rep(perm[[iter]][3],m))
  if (mean(M$classification != classtmp)<L) {
    L = mean(M$classification != classtmp)
  }
}
LR = L
MR = M

tau = as.matrix(read.csv(paste("./testdata/sim-tauhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
tau = array(tau, n)

ptau = as.matrix(read.csv(paste("./testdata/sim-ptauhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
ptau = array(ptau, c(d, n))

mu = as.matrix(read.csv(paste("./testdata/sim-muhat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
mu = array(mu, c(d, d))

sigma = as.matrix(read.csv(paste("./testdata/sim-sigmahat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
sigma = array(sigma, c(d, d, d))

pro = as.matrix(read.csv(paste("./testdata/sim-prohat",i,"-n",n,"-eps",eps,sep=""),sep = " ",header = F))
pro = array(pro, d)

write(tauinit, paste("./testdata/sim-tauinit",i,"-n",n,"-eps",eps,sep=""), ncol=1)

M$parameters$pro = pro
M$parameters$mean = t(mu)
M$parameters$variance$sigma = sigma
M$classification = tau
M$z = t(ptau)

perm = permn(1:3)
L = 3*m
for (iter in 1:length(perm)) {
  classtmp = c(rep(perm[[iter]][1],m),rep(perm[[iter]][2],m),rep(perm[[iter]][3],m))
  if (mean(tau != classtmp)<L) {
    L = mean(tau != classtmp)
  }
}
LMatlab = L


## object: output from Mclust,
## clcol: true labels, a vector of length n.
##
## ex: 
##  clcol <- rainbow(length(unique(cl))), where cl: true labels,
##  plotmclust(mout, clcol)
##
plotmclust <- function(object,clcol)
{
  require(mclust)
  
  #    object <- irisMclust
  data <- eval.parent(object$call$data)
  data <- as.matrix(data)
  G <- object$G
  p <- ncol(data)
  dimens <- seq(p)
  d <- length(dimens)
  #    cl <- as.numeric(factor(iris[,5],labels=seq(length(unique(iris[,5])))))
  
  par(mfrow=c(d,d),mar=rep(c(0.3,0.3/2),each=2),oma=c(4,4,4,4))
  for (i in seq(d)) {
    for (j in seq(d)) {
      if (i == j) {
        plot(0, 0, type = "n", xlab = "", ylab = "", axes = FALSE)
        text(0, 0, colnames(data[, dimens])[i], cex = 1.5, adj = 0.5)
        box()
      }
      else {
        coordProj(data = data, what = "classification",
                  parameters = object$parameters, z = object$z, 
                  dimens = dimens[c(j, i)], identify = FALSE, 
                  xaxt = "n", yaxt = "n",symbols=1:G)
        points(data[,c(j,i)], col=clcol,pch=object$class)
        mu <- object$parameters$mean
        sigma <- object$parameters$variance$sigma
        mu <- array(mu[dimens[c(j,i)], ], c(2, G))
        sigma <- array(sigma[dimens[c(j,i)], dimens[c(j,i)], ], c(2, 2, G))
        for (k in 1:G) mvn2plot(mu=mu[,k],sigma=sigma[,,k],k=15)
      }
      if (i == 1 && (!(j%%2))) 
        axis(3)
      if (i == d && (j%%2)) 
        axis(1)
      if (j == 1 && (!(i%%2))) 
        axis(2)
      if (j == d && (i%%2)) 
        axis(4)
    }
  }
}

clcol <- rainbow(length(unique(classesi)))[classesi]
plotmclust(MR, clcol)

plotmclust(M, clcol)
