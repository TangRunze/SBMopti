function [tauHat, muHat] = clusterX(X, nBlock, isGMM)

if (isGMM == 1)
    % Using GMM
    try
        gm = fitgmdist(X, nBlock, 'Replicates', 10);
    catch exception
        gm = fitgmdist(X, nBlock, 'Replicates', 10, 'Regularize', 0.01);
    end
    tauHat = cluster(gm, X)';
    muHat = gm.mu;
else
    % Using K-Means
    [tauHat, muHat] = kmeans(X, nBlock, 'Replicates', 10);
    tauHat = tauHat';
end