function [tauHat, muHat] = clusterX(X, nBlock, isGMM)

replicates = 5;

if (isGMM == 1)
    % Using GMM
    try
        gm = fitgmdist(X, nBlock, 'Replicates', replicates);
    catch exception
        gm = fitgmdist(X, nBlock, 'Replicates', replicates, 'Regularize', 0.01);
    end
    tauHat = cluster(gm, X)';
    muHat = gm.mu;
else
    % Using K-Means
    [tauHat, muHat] = kmeans(X, nBlock, 'Replicates', replicates);
    tauHat = tauHat';
end