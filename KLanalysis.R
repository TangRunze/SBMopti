#setwd('/Users/Runze/Dropbox/GitHub/SBMopti_KL/')
setwd('E:/Dropbox/GitHub/SBMopti_KL/')

require(ks)
require(MCMCpack)

K <- 3
n=150
diag=0.5
offdiag=0.1

nuStar <- matrix(c(0.0681,0.0681,0.7005,0.7005,0.0681,0.0681,0.0681,0.7005,0.0681),ncol=3)
# nuStar <- matrix(c(0.2453,0.6925,0.2453,0.6925,0.2453,0.2453,0.2453,0.2453,0.6925),ncol=3)
rho <- rep(1/K,K)

rVec = c(1,3,5,6,7,8,9,10,15,20)
lambdaVec = c(0.001,0.01,0.1,0.9,0.99,0.999)

kl <- array(0, c(length(rVec), length(lambdaVec)))

for (iR in 1:length(rVec)) {
	r <- rVec[iR]
	for (iLambda in 1:length(lambdaVec)) {
		lambda = lambdaVec[iLambda]
		klTmp <- c();
		for (iGraph in 1:500) {
			if ((file.exists(paste("results/KLdata/xBest-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep="")) == TRUE) & (file.exists(paste("results/KLresult/KL-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep="")) == FALSE)) {
				print(iGraph)
				xBest = as.matrix(read.csv(paste("results/KLdata/xBest-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep=""), header=F))
				kdeResult <- kde(xBest, xmin = rep(0,1,n), xmax = rep(1,1,n))
				pEmpirical <- kdeResult$estimate
				x <- kdeResult$eval.points
				klTmp <- c(klTmp,KLCalculator(x, pEmpirical, nuStar, rho, r))
				write(klTmp, file = paste("results/KLresult/KL-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep=""))
			} else if ((file.exists(paste("results/KLdata/xBest-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep="")) == TRUE) & (file.exists(paste("results/KLresult/KL-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep="")) == TRUE)) {
			  tmp <- as.matrix(read.delim(paste("results/KLresult/KL-n",n,"-diag",diag,"-offdiag",offdiag,"-r",r,"-lambda",lambda,"-adjMatrix",iGraph,sep=""), header = F, sep = " "))
			  tmp <- as.vector(t(tmp))
        nv <- !is.na(tmp)
        tmp <- tmp[nv]
        klTmp <- c(klTmp, tmp[length(tmp)])
			}
		}
		kl[iR,iLambda] <- median(klTmp)
	}
}






KLCalculator <- function(x, pEmpirical, nu, rho, r) {
	kl <- 0
	
	K <- length(x)
	
	nuTmp <- matrix(rep(0, K*(K+1)), ncol=K+1)
	nuTmp[,1:K] <- nu
	for (i in 1:K) {
		nuTmp[i,K+1] = 1 - sum(nuTmp[i,])
	}
	nu <- nuTmp
	
	nGrid <- length(x[[1]])
	area <- 1/(nGrid - 1)^K
	for (i1 in 2:(nGrid-1)) {
		x1 <- x[[1]][i1]
		for (i2 in 1:(nGrid - i1)) {
			if (i2 > 1) {
				x2 <- x[[2]][i2]
				for (i3 in 1:(nGrid - i1 - i2 + 1)) {
					if (i3 > 1) {
						x3 <- x[[3]][i3]
						pDir = rho[1]*ddirichlet(c(x1,x2,x3,1-x1-x2-x3),r*nu[1,])+rho[2]*ddirichlet(c(x1,x2,x3,1-x1-x2-x3),r*nu[2,])+rho[3]*ddirichlet(c(x1,x2,x3,1-x1-x2-x3),r*nu[3,])
						#if ((pEmpirical[i1,i2,i3] > 1e-12) & (pDir > 1e-12)) {
						if ((pEmpirical[i1,i2,i3] > 1e-12)) {
							kl <- kl + area*pEmpirical[i1,i2,i3]*log(pEmpirical[i1,i2,i3]/pDir)
						}
					}
				}
			}
		}
	}
	
	return(kl)
}

