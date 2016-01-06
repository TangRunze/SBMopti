function [tauHat, muHat] = clusterX(X, nBlock, isGMM)

replicatesGMM = 1;
replicatesKMeans = 1;

if (isGMM == 1)
    % Using GMM
    % gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM);
    % Using GMM++
    gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM, 'Start', 'plus');
    tauHat = cluster(gm, X)';
    muHat = gm.mu;
else
    % Using K-Means++
    [tauHat, muHat] = kmeans(X, nBlock, 'Replicates', replicatesKMeans);
    tauHat = tauHat';
end