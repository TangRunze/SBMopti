function [] = sbmoptisim(nVertex, nBlock, muB, epsilonInB, r, ...
    gStart, gEnd, nCore, maxIter, tol)

%% --- Quick Setting ---
% nVertex = 15;
% nBlock = 3;
% dimLatentPosition = nBlock;
% epsilonInB = 0.2;
% rho = repmat(1/nBlock, 1, nBlock);
% tol = 1e-4;
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
if (nargin < 10)
    tol = 1e-4;
end

if (nargin < 9)
    maxIter = 100;
end

if (nargin < 8)
    nCore = 1;
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
B = [0.42, 0.42; 0.42, 0.5];
rho = [0.6, 0.4];

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

%% --- Optimization Setting ---
options = optimoptions('fmincon', 'TolX', 1e-6, ...
    'MaxIter', 5000, 'MaxFunEvals', 5000, ...
    'Algorithm', 'interior-point'); %, 'GradObj', 'on');
% Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);

%% --- Parallel Computing ---
delete(gcp('nocreate'))
parpool(nCore);

% lambdaVec = [0.001 0.01 0.5 0.8 0.9 0.95 0.99 0.995 0.999];
lambdaVec = [0.001 0.01 0.1 0.9 0.99 0.999];

parfor iGraph = gStart:gEnd
    
    rng shuffle;
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    % if r = -1 then follows a SBM exactly, otherwise follows a Dirichlet
    % prior.
    [adjMatrix, nuHat0, ~, tauHat0, ~] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, ...
        rho, tauStar, r, iGraph, projectoptions);
    xHat0 = asge(adjMatrix, dimLatentPosition);
    
    % Pre-projection
    rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
        dimLatentPosition, xHat0), ...
        reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
        [], [], [], [], - ones(dimLatentPosition^2, 1), ...
        ones(dimLatentPosition^2, 1), @(x) ...
        projectconditionfun(x, nVertex, dimLatentPosition, xHat0), ...
        projectoptions);
    rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
        dimLatentPosition);
    
    % Rotate the latent positions
    xHat0 = xHat0*rotationMatrix;
    nuHat0 = nuHat0*rotationMatrix;
    
    %% --- ASGE ---
    errorRateASGE = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);
    
    for lambda = lambdaVec
        
        saveFile = ['./results/results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-r' num2str(r) '-lambda' num2str(lambda) ...
            '-adjMatrix' num2str(iGraph) '.mat'];
        
        if exist(saveFile, 'file') == 0
            %% --- Solve Optimization Problem ---
            hasConverge = 0;
            iter = 0;
            fValueBest = 0;
            xHatTmp = xHat0;
            nuHat = nuHat0;
            tauHat = tauHat0;
            
            fValueOld = objectivefun_std(reshape(xHatTmp, 1, ...
                nVertex*dimLatentPosition), adjMatrix, nuHat, lambda, ...
                nVertex, dimLatentPosition, tauHat);
            
            while (~hasConverge) && (iter < maxIter)
                iter = iter + 1;
                
                % Find the best xHatTmp based on current nuHat, tauHat.
                xHatTmp = fmincon(@(x) objectivefun_std(x, adjMatrix, ...
                    nuHat, lambda, nVertex, dimLatentPosition, tauHat), ...
                    reshape(xHatTmp, 1, nVertex*dimLatentPosition), [], ...
                    [], [], [], zeros(dimLatentPosition*nVertex, 1), ...
                    ones(dimLatentPosition*nVertex, 1), @(x) ...
                    conditionfun(x, nVertex, dimLatentPosition), options);
                xHatTmp = reshape(xHatTmp, nVertex, dimLatentPosition);
                
                % Find the best nuHat, tauHat based on current xHatTmp.
                [tauHat, nuHat] = kmeans(xHatTmp, nBlock);
                tauHat = tauHat';
                
                fValueNew = objectivefun_std(reshape(xHatTmp, 1, ...
                    nVertex*dimLatentPosition), adjMatrix, nuHat, ...
                    lambda, nVertex, dimLatentPosition, tauHat);
                
                if (fValueNew < fValueBest) || (iter == 1)
                    fValueBest = fValueNew;
                    tauBest = tauHat;
                    xBest = xHatTmp;
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
            if (hasConverge == 1)
                errorRateOpti = errorratecalculator(tauStar, tauHat, ...
                    nVertex, nBlock);
                errorRateOptiBest = errorratecalculator(tauStar, ...
                    tauBest, nVertex, nBlock);
                % Calculate the log-likelihood of the solution
                loglik = 0;
                tmpB = nuBest*nuBest';
                for i = 1:(nVertex - 1)
                    for j = (i+1):nVertex
                        loglik = loglik + adjMatrix(i, j)*...
                            log(tmpB(tauBest(i), tauBest(j))) + ...
                            (1 - adjMatrix(i, j))*...
                            log(1 - tmpB(tauBest(i), tauBest(j)));
                    end
                end
                parsave(saveFile, errorRateASGE, errorRateOpti, ...
                    errorRateOptiBest, fValueNew, fValueBest, ...
                    xBest, nuBest, tauBest, loglik);
            end
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

