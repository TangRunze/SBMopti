function [] = sbmopti(nBlock, gStart, gEnd, nCore, lambda, maxIter, tol)

%% --- Parameter Setting ---
% nVertex selects the number of vertices in the graph.
% nVertex = 150;

% nBlock selects the number of blocks in the stochastic blockmodel.
% nBlock = 3;

% dimLatentPosition selects the dimension of latent positions.
dimLatentPosition = nBlock;

% true block proportion
% rho = repmat(1/nBlock, 1, nBlock);

% epsilonInB controls the true model. The probability matrix
%       B = (0.5 - epsilonInB)*J + 2*epsilonInB*I
% epsilonInB should be inside [0, 0.5].
% epsilonInB = 0.1;

%% --- Default Parameter Setting ---
if (nargin < 7)
    tol = 1e-4;
end

if (nargin < 6)
    maxIter = 100;
end

if (nargin < 5)
    lambda = 1;
end

if (nargin < 4)
    nCore = 1;
end

if (nargin < 3)
    error('Not enough input!')
end

if ((ceil(nBlock) ~= floor(nBlock)) || (nBlock <= 0))
    error('Number of blocks should be a positive integer!')
end

if ((ceil(gStart) ~= floor(gStart)) || (ceil(gEnd) ~= floor(gEnd)) || ...
        (gStart <= 0) || (gEnd <= 0))
    error('gStart/gEnd should be positive integers!')
end

if (gStart > gEnd)
    error('gStart should be less or equal to gEnd!')
end

%% --- Parallel Computing ---
% if isempty(gcp('nocreate'))
%     parpool(nCore);
% end
delete(gcp('nocreate'))
parpool(nCore);

parfor iGraph = gStart:gEnd
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [nVertex, adjMatrix, nuHat, ~, tauHat, ~, tauStar] = ...
        datareader(nBlock, dimLatentPosition, iGraph);
    xHat = asge(adjMatrix, dimLatentPosition);
    
    %% --- ASGE ---
    errorRateASGE = errorratecalculator(tauStar, tauHat, nVertex, nBlock);
    
    saveFile = ['./results/results-SBMopti-real-graph' num2str(iGraph) ...
        '.mat'];
    if exist(saveFile, 'file') == 0
        %% --- Solve Optimization Problem ---
        options = optimoptions('fmincon', 'TolX', 1e-6, ...
            'MaxIter', 10000, 'MaxFunEvals', 10000);
        projectoptions = optimoptions('fmincon', 'TolX', 1e-6, ...
            'MaxIter', 10000, 'MaxFunEvals', 10000);
        
        hasConverge = 0;
        iter = 0;
        fValueNew = 0;
        xHatTmp = xHat;
        
        % Pre-projection
        rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
            dimLatentPosition, xHatTmp), ...
            reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
            [], [], [], [], - ones(dimLatentPosition^2, 1), ...
            ones(dimLatentPosition^2, 1), @(x) ...
            projectconditionfun(x, nVertex, dimLatentPosition, xHatTmp),...
            projectoptions);
        rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
            dimLatentPosition);
        
        % Rotation
        xHatTmp = xHatTmp*rotationMatrix;
        nuHat = nuHat*rotationMatrix;
        
        while (~hasConverge) && (iter < maxIter)
            iter = iter + 1
            fValueOld = fValueNew;
            
%             objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
%                 adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat)
            
            % Find the best xHatTmp based on current nuHat, tauHat.
            [xHatTmp, fValueNew] = fmincon(@(x) objectivefun(x, ...
                adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, ...
                tauHat), reshape(xHatTmp, 1, nVertex*dimLatentPosition),...
                [], [], [], [], zeros(dimLatentPosition*nVertex, 1), ...
                ones(dimLatentPosition*nVertex, 1), @(x) ...
                conditionfun(x, nVertex, dimLatentPosition), options);
            xHatTmp = reshape(xHatTmp, nVertex, dimLatentPosition);
            
            % Convergence Checking
            if (abs(fValueNew - fValueOld) < tol*fValueOld) && (iter > 1)
                hasConverge = 1;
            end
            
            [tauHat, nuHat] = kmeans(xHatTmp, nBlock);
            tauHat = tauHat';
        end
        
        errorRateOpti = errorratecalculator(tauStar, tauHat, nVertex, ...
            nBlock);
        
        %% --- Save Results ---
        if (hasConverge == 1)
            parsave(saveFile, errorRateASGE, errorRateOpti);
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

