function [tauHat, muHat] = clusterX(X, nBlock, isGMM)

replicatesGMM = 1;
replicatesKMeans = 1;

if (isGMM == 1)
    % Using GMM
    % gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM);
    
    % Using GMM++
    % gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM, 'Start', 'plus');
    
    % Covariance matrices = I
    % gm = fitgmdist(X, nBlock, 'SharedCovariance', true);
    
    % Initialized using K-Means++
    [~, muHat] = kmeans(X, nBlock, 'Replicates', replicatesKMeans);
    S.mu = muHat;
    S.Sigma = ones(1, size(X, 2));
    S.ComponentProportion = repmat(1/nBlock, 1, nBlock);
    gm = fitgmdist(X, nBlock, 'CovarianceType', 'diagonal', ...
        'SharedCovariance', true, 'Start', S);
    
    tauHat = cluster(gm, X)';
    muHat = gm.mu;
else
    % Using K-Means++
    [tauHat, muHat] = kmeans(X, nBlock, 'Replicates', replicatesKMeans);
    tauHat = tauHat';
end