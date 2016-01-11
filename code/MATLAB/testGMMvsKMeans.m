close all
clear all


%% Load data
% load('data/testGMMvsKMeans.mat');

r = -1;

nBlock = 2;
dimLatentPosition = 2;
ond = 0.3;
ofd = 0.2;

nVertex = 150;
rho = ones(nBlock,1)/nBlock;
nB = rho(1)*nVertex;
tauStar = [ones(nB,1)', 2*ones(nB,1)']; %3*ones(nB,1)]';
% B=[ond ofd ofd; ofd ond ofd; ofd ofd ond];
B = [ond ofd; ofd ond];

% projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter',100, ...
%     'MaxFunEvals', 5000);

%% sample graphs and estimate stuff

nmc = 100;
errorRateASGE_k = nan(nmc,1);
errorRateASGE_g = nan(nmc,1);

for i=1:nmc
    
    [~, adjMatrixDA] =  datagenerator(nVertex, nBlock, ...
        dimLatentPosition, B, rho, tauStar, r, i);
    
    % Calculation
    % Calculate estimated latent positions
    xHat = asge(adjMatrixDA, dimLatentPosition);
    
    % Cluster using K-Means
    [tauHat_k, ~] = clusterX(xHat, nBlock, 0);
    errorRateASGE_k(i) = errorratecalculator(tauStar, tauHat_k, ...
        nVertex, nBlock);
    
    % Cluster using GMM
    [tauHat_g, ~] = clusterX(xHat, nBlock, 1);
    errorRateASGE_g(i) = errorratecalculator(tauStar, tauHat_g, ...
        nVertex, nBlock);
    
end
i = nmc;

hist(errorRateASGE_g - errorRateASGE_k)
title('histogram of error rate of GMM - error rate of K-Means');
xlabel('error rate of GMM - error rate of K-Means');
ylabel('Count over 100 replicates');
[mean(errorRateASGE_k) mean(errorRateASGE_g)]

%%
eps = 10^-4;
figure(1), clf, hold all
plot(errorRateASGE_g + eps,errorRateASGE_k + eps, '.', 'markersize', 12)
xlabel('kmeans'), ylabel('gmm')
m = min(min(errorRateASGE_k), min(errorRateASGE_g));
plot([m, 0.6], [m, 0.6])
set(gca, 'xscale', 'log', 'yscale', 'log')
%%
figure(2), clf, hold all
plot(errorRateASGE_g + eps), set(gca, 'yscale', 'log')
plot(errorRateASGE_k + eps), set(gca, 'yscale', 'log')
legend('gmm', 'kmeans')

figure(3), clf, 
hist(errorRateASGE_g - errorRateASGE_k, 50)

%%

i = i + 1;
[~, adjMatrixDA, nuStar] =  datagenerator(nVertex, nBlock, ...
    dimLatentPosition, B, rho, tauStar, r, i);

% Calculation
% Calculate estimated latent positions
xHat = asge(adjMatrixDA, dimLatentPosition);

% Cluster using K-Means
[tauHat_k, nuHat_k] = clusterX(xHat, nBlock, 0);
errorRateASGE_k(i) = errorratecalculator(tauStar, tauHat_k, nVertex, nBlock);

% Cluster using GMM
[tauHat_g, nuHat_g] = clusterX(xHat, nBlock, 1);
errorRateASGE_g(i) = errorratecalculator(tauStar, tauHat_g, nVertex, nBlock);


[errorRateASGE_k(i), errorRateASGE_g(i)]

figure(4), clf, plot2Dxhat(xHat, nuHat_k, tauStar, nuStar, nuHat_g); 
title(['kmeans=', num2str(errorRateASGE_k(i)), ' gmm=', ...
    num2str(errorRateASGE_g(i))])








%% pdf

i = i+1;
[adjMatrix, adjMatrixDA, nuStar] =  datagenerator(nVertex, nBlock, ...
    dimLatentPosition, B, rho, tauStar, r, i);

% Calculation
% Calculate estimated latent positions
xHat = asge(adjMatrixDA, dimLatentPosition);

% Cluster using K-Means
[tauHat_k, nuHat_k] = clusterX(xHat, nBlock, 0);
errorRateASGE_k(i) = errorratecalculator(tauStar, tauHat_k, nVertex, nBlock);

% Cluster using GMM
% [tauHat_g, nuHat_g] = clusterX(xHat, nBlock, 1);
% errorRateASGE_g(i) = errorratecalculator(tauStar, tauHat_g, nVertex, nBlock);
gm = fitgmdist(xHat, nBlock, 'Start', 'plus');
tauHat_g = cluster(gm, xHat)';
nuHat_g = gm.mu;
sigmaHat_g = gm.Sigma;
proportionHat_g = gm.ComponentProportion;

errorRateASGE_g(i) = errorratecalculator(tauStar, tauHat_g, nVertex, nBlock);


[errorRateASGE_k(i), errorRateASGE_g(i)]

figure(5), clf, plot3DxhatPDF(xHat, tauStar, nuStar, nuHat_k, nuHat_g, ...
    sigmaHat_g, proportionHat_g); 
title(['kmeans=', num2str(errorRateASGE_k(i)), ' gmm=', ...
    num2str(errorRateASGE_g(i))])



% Check GMM clustering
figure(6), clf, plot3DxhatPDF(xHat, tauHat_g, nuStar, nuHat_k, nuHat_g, ...
    sigmaHat_g, proportionHat_g); 



% Check GMM clustering
figure(7), clf, plot3DxhatPDF(xHat, tauHat_k, nuStar, nuHat_k, nuHat_g, ...
    sigmaHat_g, proportionHat_g); 



figure(8), clf, plot2Dxhatheatmap(xHat, nuHat_k, tauStar, nuStar, nuHat_g); 
title(['kmeans=', num2str(errorRateASGE_k(i)), ' gmm=', ...
    num2str(errorRateASGE_g(i))])

