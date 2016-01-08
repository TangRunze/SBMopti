close all
clear all

%% --- Simulation ---

r = 20;
iGraph = 2;

% pathExtra = '/';
pathExtra = '/latest/';

% Parameter Setting
nVertex = 150;

rho = [1/3, 1/3, 1/3];

B = [0.4, 0.2, 0.2; 0.2, 0.4, 0.2; 0.2, 0.2, 0.4];

nBlock = 3;
dimLatentPosition = 3;
isGMM = 0;

tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    100, 'MaxFunEvals', 5000);




[adjMatrix, adjMatrixDA] = ...
    datagenerator(nVertex, nBlock, dimLatentPosition, B, ...
    rho, tauStar, r, iGraph, projectoptions);

% Calculate estimated latent positions
xHat0 = asge(adjMatrixDA, dimLatentPosition);

% Cluster using K-Means
[tauHat0, nuHat0] = clusterX(xHat0, nBlock, 0);
errorRateASGE0 = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);

% Cluster using GMM
[tauHat1, nuHat1] = clusterX(xHat0, nBlock, 1);
errorRateASGE1 = errorratecalculator(tauStar, tauHat1, nVertex, nBlock);

[errorRateASGE0, errorRateASGE1]

% plot3Dxhat(xHat0, nuHat0, tauStar);
% plot3Dxhat(xHat0, nuHat1, tauStar);










