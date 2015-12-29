function [] = sbmoptisim(nVertex, nBlock, muB, epsilonInB, r, ...
    gStart, gEnd, isGMM, nCore, maxIter, tol)

%% --- Quick Setting ---
% nVertex = 15;
% nBlock = 3;
% dimLatentPosition = nBlock;
% epsilonInB = 0.2;
% rho = repmat(1/nBlock, 1, nBlock);
% tol = 1e-4;sbmoptisim(150, 3, 0.3, 0.1, 0, 1, 200, 0, 12)
% maxIter = 100;
% muB = 0.3;
% r = 100;
% B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
% tauStar = [];
% nVectorStar = nVertex*rho;
% for i = 1:nBlock
%     tauStar = [tauStar, i*ones(1, nVectorStar(i))];
% end
% options = optimoptions('fmincon', 'TolX', 1e-6, ...
%     'MaxIter', 10000, 'MaxFunEvals', 10000, ...
%     'Algorithm', 'interior-point');
% projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
%     10000, 'MaxFunEvals', 10000);

%% --- Parameter Setting ---
% nVertex selects the number of vertices in the graph.
% nVertex = 150;

% nBlock selects the number of blocks in the stochastic blockmodel.
% nBlock = 3;

% dimLatentPosition selects the dimension of latent positions.
dimLatentPosition = nBlock;

% true block proportion
rho = repmat(1/nBlock, 1, nBlock);

%% --- Default Parameter Setting ---
if (nargin < 11)
    tol = 1e-4;
end

if (nargin < 10)
    maxIter = 100;
end

if (nargin < 9)
    nCore = 1;
end

if (nargin < 8)
    isGMM = 0;
end

if (nargin < 7)
    error('Not enough input!')
end

if ((ceil(nVertex) ~= floor(nVertex)) || (nVertex <= 0))
    error('Number of vertices should be a positive integer!')
end

if ((ceil(nBlock) ~= floor(nBlock)) || (nBlock <= 0))
    error('Number of blocks should be a positive integer!')
end

if ((muB - epsilonInB < 0) || (muB + epsilonInB > 1))
    error('Probability matrix invalid!')
end

if ((ceil(gStart) ~= floor(gStart)) || (ceil(gEnd) ~= floor(gEnd)) || ...
        (gStart <= 0) || (gEnd <= 0))
    error('gStart/gEnd should be positive integers!')
end

if (gStart > gEnd)
    error('gStart should be less or equal to gEnd!')
end

% block probability matrix
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

%%% --- For Minh's Model --- %%%
% B = [0.42, 0.42; 0.42, 0.5];
% rho = [0.6, 0.4];

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

%% --- Optimization Setting ---
options = optimoptions('fmincon', 'TolX', 1e-4, ...
    'MaxIter', 100, 'MaxFunEvals', 5000, ...
    'Algorithm', 'interior-point'); %, 'GradObj', 'on');
% Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    100, 'MaxFunEvals', 5000);

%% --- Parallel Computing ---
delete(gcp('nocreate'))
parpool(nCore);

% lambdaVec = [0.001 0.01 0.5 0.8 0.9 0.95 0.99 0.995 0.999];
lambdaVec = [0.001 0.01 0.1 0.5 0.9 0.99 0.999];

parfor iGraph = gStart:gEnd
    
    rng shuffle;
    
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
    
    % Get a better initialization (most positive)
    rotationMatrix0 = eye(dimLatentPosition);
    for iD = 1:dimLatentPosition
        if sum(xHat0(:, iD) > 0) < nVertex/2
            rotationMatrix0(iD, iD) = -1;
        end
    end
    xHat0 = xHat0*rotationMatrix0;
    nuHat0 = nuHat0*rotationMatrix0;
    
    % Optimize over mean
    xHat0mean = mean(xHat0);
    rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
        dimLatentPosition, xHat0mean), ...
        reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
        [], [], [], [], - ones(dimLatentPosition^2, 1), ...
        ones(dimLatentPosition^2, 1), @(x) ...
        projectconditionfun(x, 1, dimLatentPosition, xHat0mean),...
        projectoptions);
    rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
        dimLatentPosition);
    
    % Rotation
    xHat0 = xHat0*rotationMatrix;
    nuHat0 = nuHat0*rotationMatrix;
    
    % Make nuHat0 feasible
    nuHat0(nuHat0 < 0) = 1e-6;
    rowSumTmp = sqrt(sum(nuHat0.*nuHat0, 2));
    nuHat0(rowSumTmp > 1, :) = nuHat0(rowSumTmp > 1, :)./ ...
        repmat(rowSumTmp(rowSumTmp > 1), 1, dimLatentPosition);

    % Make xHat0 feasible
    xHat0(xHat0 < 0) = 1e-6;
    rowSumTmp = sqrt(sum(xHat0.*xHat0, 2));
    xHat0(rowSumTmp > 1, :) = xHat0(rowSumTmp > 1, :)./ ...
        repmat(rowSumTmp(rowSumTmp > 1), 1, dimLatentPosition);

    
    for lambda = lambdaVec
        
        saveFile = ['./results/results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-GMM' num2str(isGMM) '-r' num2str(r) ...
            '-lambda' num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'];
        
        if exist(saveFile, 'file') == 0
            %% --- Solve Optimization Problem ---
            hasConverge = 0;
            iter = 0;
            xHat = xHat0;
            nuHat = nuHat0;
            tauHat = tauHat0;
            
            fValueOld = objectivefun_std(reshape(xHat, 1, ...
                nVertex*dimLatentPosition), adjMatrix, nuHat, lambda, ...
                nVertex, dimLatentPosition, tauHat);
            
            fValueBest = fValueOld;
            tauBest = tauHat;
            xBest = xHat;
            nuBest = nuHat;
            
            while (~hasConverge) && (iter < maxIter)
                iter = iter + 1;
                
                % Find the best xHat based on current nuHat, tauHat.
                xHat = fmincon(@(x) objectivefun_std(x, adjMatrix, ...
                    nuHat, lambda, nVertex, dimLatentPosition, tauHat), ...
                    reshape(xHat, 1, nVertex*dimLatentPosition), [], ...
                    [], [], [], zeros(dimLatentPosition*nVertex, 1), ...
                    ones(dimLatentPosition*nVertex, 1), @(x) ...
                    conditionfun(x, nVertex, dimLatentPosition), options);
                xHat = reshape(xHat, nVertex, dimLatentPosition);
                
                % Find the best nuHat, tauHat based on current xHat.
                [tauHat, nuHat] = clusterX(xHat, nBlock, isGMM);
                
                fValueNew = objectivefun_std(reshape(xHat, 1, ...
                    nVertex*dimLatentPosition), adjMatrix, nuHat, ...
                    lambda, nVertex, dimLatentPosition, tauHat);
                
                if (fValueNew < fValueBest)
                    fValueBest = fValueNew;
                    tauBest = tauHat;
                    xBest = xHat;
                    nuBest = nuHat;
                end
                
                % Convergence Checking
                if (abs(fValueNew - fValueOld) < tol*fValueOld) && ...
                        (iter > 1)
                    hasConverge = 1;
                end
                
                fValueOld = fValueNew;
            end
            
            %% --- Save Results ---
            if (hasConverge == 1) || (iter == maxIter)
                % Calculate the error rate
                errorRateOpti = errorratecalculator(tauStar, tauBest, ...
                    nVertex, nBlock);
                % Calculate the log-likelihood of the solution
                loglikASGE = loglikcalculator(adjMatrix, tauHat0, nBlock);
                loglikOpti = loglikcalculator(adjMatrix, tauBest, nBlock);
                
                parsavesim(saveFile, errorRateASGE, errorRateOpti, xBest, ...
                    nuBest, tauBest, xHat0, nuHat0, tauHat0, loglikASGE,...
                    loglikOpti);
            end
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

