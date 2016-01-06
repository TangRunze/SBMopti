close all
clear all


%% Load data
load('testGMMvsKMeans.mat');

%% Calculation

% Calculate estimated latent positions
xHat0 = asge(adjMatrixDA, dimLatentPosition);

% Cluster using K-Means
[tauHat0, nuHat0] = clusterX(xHat0, nBlock, 0);
errorRateASGE0 = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);

% Cluster using GMM
[tauHat1, nuHat1] = clusterX(xHat0, nBlock, 1);
errorRateASGE1 = errorratecalculator(tauStar, tauHat1, nVertex, nBlock);

[errorRateASGE0, errorRateASGE1]


%% Plot
% plot3Dxhat(xHat0, nuHat0, tauStar);
% plot3Dxhat(xHat0, nuHat1, tauStar);
