
clear all

%% --- Parameter Setting ---
% nVertex selects the number of vertices in the graph.
nVertex = 150;

% nBlock selects the number of blocks in the stochastic blockmodel.
nBlock = 3;

% dimLatentPosition selects the dimension of latent positions.
dimLatentPosition = nBlock;

% true block proportion
rho = repmat(1/nBlock, 1, nBlock);

% epsilonInB controls the true model. The probability matrix
%       B = (0.5 - epsilonInB)*J + 2*epsilonInB*I
% epsilonInB should be inside [0, 0.5].
epsilonInB = 0.1;

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

%% --- Generate/Read Data ---
% Generate data if there does not exist one, otherwise read the
% existing data.
[adjMatrix, nuHat, sigmaHat, tauHat, pTauHat] = ...
    datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, ...
    epsilonInB, 1);
xHat = asge(adjMatrix, dimLatentPosition);

%% --- ASGE ---
errorRateASGE = errorratecalculator(tauStar, tauHat, nVertex, nBlock);

%% --- Solve Optimization Problem ---

lambda = 1;
tol = 1e-4;
maxIter = 100;
options = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', 10000, ...
    'MaxFunEvals', 10000);

projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', 100000,...
    'MaxFunEvals', 100000);

hasConverge = 0;
iter = 0;
fValueNew = 0;
xHatTmp = xHat;

% Pre-projection
rotationMatrix = fmincon(@(x) projectobjectivefun(x, dimLatentPosition, ...
    xHatTmp), reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
    [], [], [], [], - ones(dimLatentPosition^2, 1), ...
    ones(dimLatentPosition^2, 1), @(x) ...
    projectconditionfun(x, nVertex, dimLatentPosition, xHatTmp), ...
    projectoptions);
rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
    dimLatentPosition);

% Rotation
xHatTmp = xHatTmp*rotationMatrix;
nuHat = nuHat*rotationMatrix;

while (~hasConverge) && (iter < maxIter)
    iter = iter + 1
    fValueOld = fValueNew;
    
    objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
        adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat)
    
    % Find the best xHatTmp based on current nuHat, tauHat.
%     [xHatTmp, fValueNew] = fmincon(@(x) objectivefun(x, adjMatrix, ...
%         nuHat, lambda, nVertex, dimLatentPosition, tauHat), ...
%         reshape(xHatTmp, 1, nVertex*dimLatentPosition), [], [], [], [], ...
%         - ones(dimLatentPosition*nVertex, 1), ...
%         ones(dimLatentPosition*nVertex, 1), @(x) ...
%         conditionfun(x, nVertex, dimLatentPosition), options);
    [xHatTmp, fValueNew] = fmincon(@(x) objectivefun(x, adjMatrix, ...
        nuHat, lambda, nVertex, dimLatentPosition, tauHat), ...
        reshape(xHatTmp, 1, nVertex*dimLatentPosition), [], [], [], [], ...
        zeros(dimLatentPosition*nVertex, 1), ...
        ones(dimLatentPosition*nVertex, 1), @(x) ...
        conditionfun(x, nVertex, dimLatentPosition), options);
    xHatTmp = reshape(xHatTmp, nVertex, dimLatentPosition);
    
    fValueNew
    
%     objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
%         adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat)
    
    % Convergence Checking
    if (abs(fValueNew - fValueOld) < tol*fValueOld) && (iter > 1)
        hasConverge = 1;
    end
    
%     gm = fitgmdist(xHatTmp, nBlock, 'Replicates', 20);
%     tauHat = cluster(gm, xHatTmp)';
%     if (norm(gm.mu - nuHat) < tol)
%         hasConverge = 1;
%     end
%     nuHat = gm.mu;
    nuHatOld = nuHat;
    [tauHat, nuHat] = kmeans(xHatTmp, nBlock);
    tauHat = tauHat';
%     if (norm(nuHatOld - nuHat) < tol*norm(nuHatOld))
%         hasConverge = 1;
%     end

end

errorRateOpti = errorratecalculator(tauStar, tauHat, nVertex, nBlock);



