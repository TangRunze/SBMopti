function [] = sbmopti_k(nBlock, gStart, gEnd, isGMM, nCore)

%% --- Parameter Setting ---
% dimLatentPosition selects the dimension of latent positions.
dimLatentPosition = nBlock;

tol = 1e-4;
maxIter = 100;

%% --- Default Parameter Setting ---
if (nargin < 5)
    nCore = 1;
end

if (nargin < 4)
    isGMM = 0;
end

if (nargin < 3)
    error('Not enough input!')
end

if ((ceil(nBlock) ~= floor(nBlock)) || (nBlock <= 0))
    error('Number of blocks should be a positive integer!')
end

if ((ceil(gStart) ~= floor(gStart)) || (ceil(gEnd) ~= floor(gEnd)) || ...
        (gStart < 0) || (gEnd < 0))
    error('gStart/gEnd should be positive integers!')
end

if (gStart > gEnd)
    error('gStart should be less or equal to gEnd!')
end

%% --- Optimization Setting ---
% options = optimoptions('fmincon', 'TolX', 1e-4, ...
%     'MaxIter', 100, 'MaxFunEvals', 10000, 'Display', 'iter');
% projectoptions = optimoptions('fmincon', 'TolX', 1e-6, ...
%     'MaxIter', 100, 'MaxFunEvals', 10000, 'Display', 'iter');

options = optimoptions('fmincon', 'TolX', 1e-4, ...
    'MaxIter', 100, 'MaxFunEvals', 5000);
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, ...
    'MaxIter', 100, 'MaxFunEvals', 5000);

%% --- Parallel Computing ---
% if isempty(gcp('nocreate'))
%     parpool(nCore);
% end
delete(gcp('nocreate'))
parpool(nCore);

% lambdaVec = [0.001 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.999];
% lambdaVec = [0.001 0.01 0.1 0.5 0.9];
lambdaVec = [0.5];

for iGraph = gStart:gEnd
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [nVertex, adjMatrix, adjMatrixDA, ~, tauStar] = ...
        datareader_wiki(iGraph);
    
    xHat0 = asge(adjMatrixDA, dimLatentPosition);
    
    [tauHat0, nuHat0] = clusterX(xHat0, nBlock, isGMM);
    
    [tauHat0c, ~] = clusterX(xHat0, 3, isGMM);
    errorRateASGE = errorratecalculator(tauStar, tauHat0c, nVertex, 3);
    
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
    nuHat0(rowSumTmp > 1, :) = nuHat0(rowSumTmp > 1, :)./repmat(rowSumTmp(rowSumTmp > 1), 1, dimLatentPosition);

    % Make xHat0 feasible
    xHat0(xHat0 < 0) = 1e-6;
    rowSumTmp = sqrt(sum(xHat0.*xHat0, 2));
    xHat0(rowSumTmp > 1, :) = xHat0(rowSumTmp > 1, :)./repmat(rowSumTmp(rowSumTmp > 1), 1, dimLatentPosition);
    
    for iLambda = 1:length(lambdaVec)
        lambda = lambdaVec(iLambda);
        
        saveFile = ['./results/results-SBMopti-real-graph' ...
            num2str(iGraph) '-dim' num2str(dimLatentPosition) ...
            '-GMM' num2str(isGMM) '-lambda' num2str(lambda) '.mat'];
        
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
                [tauBestc, ~] = clusterX(xBest, 3, isGMM);
                errorRateOpti = errorratecalculator(tauStar, tauBestc, ...
                    nVertex, 3);
                
                % Calculate the log-likelihood of the solution
                loglikASGE = loglikcalculator(adjMatrix, tauHat0, nBlock);
                loglikASGEc = loglikcalculator(adjMatrix, tauHat0c, nBlock);
                loglikOpti = loglikcalculator(adjMatrix, tauBest, nBlock);
                loglikOptic = loglikcalculator(adjMatrix, tauBestc, nBlock);
                
                parsave(saveFile, errorRateASGE, errorRateOpti, xBest, ...
                    nuBest, tauBest, tauBestc, xHat0, nuHat0, tauHat0, ...
                    tauHat0c, loglikASGE, loglikASGEc, loglikOpti, ...
                    loglikOptic);
            end
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

