function [] = sbmoptisim_block(nVertex, nBlock, epsilonInB, delta, ...
    gStart, gEnd, nCore, lambda1, maxIter, tol)

%% --- Parameter Setting ---
% nVertex selects the number of vertices in the graph.
% nVertex = 150;

% nBlock selects the number of blocks in the stochastic blockmodel.
% nBlock = 3;

% dimLatentPosition selects the dimension of latent positions.
dimLatentPosition = nBlock;

% true block proportion
rho = repmat(1/nBlock, 1, nBlock);

% epsilonInB controls the true model. The probability matrix
%       B = (0.5 - epsilonInB)*J + 2*epsilonInB*I
% epsilonInB should be inside [0, 0.5].
% epsilonInB = 0.1;

%% --- Default Parameter Setting ---
if (nargin < 10)
    tol = 1e-4;
end

if (nargin < 9)
    maxIter = 100;
end

if (nargin < 8)
    lambda1 = 1;
end

if (nargin < 7)
    nCore = 1;
end

if (nargin < 6)
    error('Not enough input!')
end

if ((ceil(nVertex) ~= floor(nVertex)) || (nVertex <= 0))
    error('Number of vertices should be a positive integer!')
end

if ((ceil(nBlock) ~= floor(nBlock)) || (nBlock <= 0))
    error('Number of blocks should be a positive integer!')
end

if ((epsilonInB < 0) || (epsilonInB > 0.5))
    error('Epsilon should be inside [0, 0.5]!')
end

if ((ceil(gStart) ~= floor(gStart)) || (ceil(gEnd) ~= floor(gEnd)) || ...
        (gStart <= 0) || (gEnd <= 0))
    error('gStart/gEnd should be positive integers!')
end

if (gStart > gEnd)
    error('gStart should be less or equal to gEnd!')
end

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

%% --- Parallel Computing ---
% if isempty(gcp('nocreate'))
%     parpool(nCore);
% end
delete(gcp('nocreate'))
parpool(nCore);

% lambdaVec = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 10, 20, 50, 100];
lambdaVec = 0.1;

parfor iProbMatrix = gStart:gEnd
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [pMatrix, nuHat0, ~, tauHat0, ~] = ...
        datagenerator_block(nVertex, nBlock, dimLatentPosition, B, rho, ...
        tauStar, epsilonInB, delta, iProbMatrix);
    xHat = asge(pMatrix, dimLatentPosition);
    
    %% --- ASGE ---
    errorRateASGE = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);
    
    for lambda = lambdaVec
        
        saveFile = ['./results/results-SBMopti-Block-sim-n' num2str(nVertex)...
            '-eps' num2str(epsilonInB) '-delta' num2str(delta) ...
            '-lambda' num2str(lambda) '-pmatrix' num2str(iProbMatrix) '.mat'];
        
        if exist(saveFile, 'file') == 0
            
            % scatterdata = zeros(3,2);
            
            %% --- Solve Optimization Problem ---
            options = optimoptions('fmincon', 'TolX', 1e-6, ...
                'MaxIter', 10000, 'MaxFunEvals', 10000, ...
                'Algorithm', 'interior-point', 'GradObj', 'on');
            % Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
            projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
                10000, 'MaxFunEvals', 10000);
            
            hasConverge = 0;
            iter = 0;
            fValueBest = 0;
            xHatTmp = xHat;
            nuHat = nuHat0;
            tauHat = tauHat0;
            
            % Pre-projection
            rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
                dimLatentPosition, xHatTmp), ...
                reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
                [], [], [], [], - ones(dimLatentPosition^2, 1), ...
                ones(dimLatentPosition^2, 1), @(x) ...
                projectconditionfun(x, nVertex, dimLatentPosition, xHatTmp), ...
                projectoptions);
            rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
                dimLatentPosition);
            
            % Rotation
            xHatTmp = xHatTmp*rotationMatrix;
            nuHat = nuHat*rotationMatrix;
            
            fValueOld = objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
            % errorOld = errorRateASGE;
            
            while (~hasConverge) && (iter < maxIter)
                iter = iter + 1
                
                %             objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                %                 pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat)
                
                % Find the best xHatTmp based on current nuHat, tauHat.
                [xHatTmp, fValueNew] = fmincon(@(x) objectivefun(x, ...
                    pMatrix, nuHat, lambda, nVertex, dimLatentPosition, ...
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
                
                fValueTmp = objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                    pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
                if (fValueTmp < fValueBest) || (iter == 1)
                    fValueBest = fValueTmp;
                    tauBest = tauHat;
                end
                
                % Save for scatter plot
                % errorNew = errorratecalculator(tauStar, tauHat, ...
                %     nVertex, nBlock);
                % scatterdata(1, iter) = (fValueOld - fValueNew)/fValueOld;
                % scatterdata(2, iter) = (fValueOld - fValueTmp)/fValueOld;
                % scatterdata(3, iter) = errorOld - errorNew;
                % errorOld = errorNew;
                
                fValueOld = fValueNew;
            end
            
            %% --- Save Results ---
            if (hasConverge == 1)
                errorRateOpti = errorratecalculator(tauStar, tauHat, ...
                    nVertex, nBlock);
                errorRateOptiBest = errorratecalculator(tauStar, tauBest, ...
                    nVertex, nBlock);
                % parsave(saveFile, errorRateASGE, errorRateOpti, errorRateOptiBest, scatterdata);
                parsave(saveFile, errorRateASGE, errorRateOpti, ...
                    errorRateOptiBest, fValueNew, fValueBest);
            end
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

