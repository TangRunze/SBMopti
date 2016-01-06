function [tauHat, muHat] = clusterX(X, nBlock, isGMM)

replicatesGMM = 1;
replicatesKMeans = 1;

if (isGMM == 1)
    % Using GMM
%     try
%         gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM);
%     catch exception
%         gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM, 'Regularize', 1e-12);
%     end
    gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM);
%     gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM, 'Start', 'plus');
    tauHat = cluster(gm, X)';
    muHat = gm.mu;
else
    % Using K-Means
    [tauHat, muHat] = kmeans(X, nBlock, 'Replicates', replicatesKMeans);
    tauHat = tauHat';
end