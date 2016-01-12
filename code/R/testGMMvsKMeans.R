
asge = function(Ai,d) {
  # this function generates the adjacency spectral graph embedding
  S = eigen(Ai)
  Xhat = S$vectors[,1:d] %*% diag(sqrt(S$values[1:d]))
}

library(MASS);
library(stats);
library(LICORS);
library(mclust);

nVertex = 150;
nBlock = 2;
dimLatentPosition = 2;
rho = c(0.5, 0.5);
t = 1;
mu1 = c(-t, -t);
mu2 = c(t, t);
sigma1 = matrix(c(1, 0, 0, 1), ncol = 2);
sigma2 = matrix(c(1, 0, 0, 2), ncol = 2);
xHat = matrix(rep(0, nVertex*dimLatentPosition), ncol = dimLatentPosition);

nmc = 100;

errorRateASGE_k = rep(NaN, nmc);
errorRateASGE_g = rep(NaN, nmc);

for (i in 1:nmc) {
  tauStar = (runif(nVertex) > rho[1]) + 1;
  nv1 = (tauStar == 1);
  nv2 = (tauStar == 2);
  xHat[nv1, ] = mvrnorm(sum(nv1), mu1, sigma1);
  xHat[nv2, ] = mvrnorm(sum(nv2), mu2, sigma2);
  
  # Cluster using K-Means
  cl_kmeans <- kmeans(xHat, nBlock);
  tauHat_k <- cl_kmeans$cluster;
  
  errorRateASGE_k[i] = min(sum(tauHat_k != tauStar), sum(tauHat_k != (3 - tauStar)));
  
  # kmeanspp(data, k = 2, start = "random", iter.max = 100, nstart = 10, ...)
  
  # Cluster using GMM
  # cl_GMM <- Mclust(xHat, nBlock);
  cl_GMM <- Mclust(xHat, nBlock, modelNames=c("EII"));
  tauHat_g <- cl_GMM$classification;
  
  errorRateASGE_g[i] = min(sum(tauHat_g != tauStar), sum(tauHat_g != (3 - tauStar)));
  
}


hist(errorRateASGE_g - errorRateASGE_k)

c(mean(errorRateASGE_k), mean(errorRateASGE_g))
