function [tauHat, muHat] = clusterX(X, nBlock, isGMM)

replicatesGMM = 5;
replicatesKMeans = 5;

if (isGMM == 1)
    % Using GMM
    gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM);
    
    % Using GMM++
%     gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM, 'Start', 'plus');
    
    % Initialized using K-Means++, constrained
%     [~, muHat] = kmeans(X, nBlock, 'Replicates', replicatesKMeans);
%     S.mu = muHat;
%     S.Sigma = ones(1, size(X, 2));
%     S.ComponentProportion = repmat(1/nBlock, 1, nBlock);
%     gm = fitgmdist(X, nBlock, 'CovarianceType', 'diagonal', ...
%         'SharedCovariance', true, 'Start', S);
    
    % GMM++, constrained
%     gm = fitgmdist(X, nBlock, 'Replicates', replicatesGMM, ...
%         'CovarianceType', 'diagonal', 'SharedCovariance', true, ...
%         'Start', 'plus');
    
    tauHat = cluster(gm, X)';
    muHat = gm.mu;
else
    % Using K-Means++
    [tauHat, muHat] = kmeans(X, nBlock, 'Replicates', replicatesKMeans);
    tauHat = tauHat';
end