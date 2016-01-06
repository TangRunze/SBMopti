
r = 20;

isGMM = 1;

maxIter = 200;
nVertex = 150;
nBlock = 3;
B = [0.4, 0.2, 0.2; 0.2, 0.4, 0.2; 0.2, 0.2, 0.4];
rho = [1/3, 1/3, 1/3];
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

options = optimoptions('fmincon', 'TolX', 1e-4, ...
    'MaxIter', 100, 'MaxFunEvals', 5000, ...
    'Algorithm', 'interior-point'); %, 'GradObj', 'on');
% Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    100, 'MaxFunEvals', 5000);

errorASGEVec = zeros(1, maxIter);

for iGraph = 1:maxIter
    
    iGraph
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    % if r = -1 then follows a SBM exactly, otherwise follows a Dirichlet
    % prior.
    [adjMatrix, adjMatrixDA] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, ...
        rho, tauStar, r, iGraph, projectoptions);
    
    xHat0 = asge(adjMatrixDA, dimLatentPosition);
    
    [tauHat0, nuHat0] = clusterX(xHat0, nBlock, isGMM);
    errorRateASGE = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);
    
    errorASGEVec(iGraph) = errorRateASGE;
end

mean(errorASGEVec)

