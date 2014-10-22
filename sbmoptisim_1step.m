function [] = sbmoptisim_1step(nVertex, nBlock, epsilonInB, gStart, gEnd, ...
    nCore, lambda, maxIter, tol)

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
if (nargin < 9)
    tol = 1e-4;
end

if (nargin < 8)
    maxIter = 100;
end

if (nargin < 7)
    lambda = 1;
end

if (nargin < 6)
    nCore = 1;
end

if (nargin < 5)
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

parfor iGraph = gStart:gEnd
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [adjMatrix, nuHat, ~, tauHat, ~] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, ...
        epsilonInB, iGraph);
    xHat = asge(adjMatrix, dimLatentPosition);
    
    %% --- ASGE ---
    errorRateASGE = errorratecalculator(tauStar, tauHat, nVertex, nBlock);
    
    saveFile = ['./results/results-SBMopti-sim-n' num2str(nVertex) ...
        '-eps' num2str(epsilonInB) '-graph' num2str(iGraph) '.mat'];
    if exist(saveFile, 'file') == 0
        
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
        
        %% --- Pre-projection Method 1 ---
        rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
            dimLatentPosition, xHatTmp), ...
            reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
            [], [], [], [], - ones(dimLatentPosition^2, 1), ...
            ones(dimLatentPosition^2, 1), @(x) ...
            projectconditionfun(x, nVertex, dimLatentPosition, xHatTmp), ...
            projectoptions);
        rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
            dimLatentPosition);
        
        rotationMatrix = fmincon(@(x) 0, ...
            reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
            [], [], [], [], - ones(dimLatentPosition^2, 1), ...
            ones(dimLatentPosition^2, 1), @(x) ...
            projectconditionfun(x, nVertex, dimLatentPosition, xHatTmp), ...
            projectoptions);
        rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
            dimLatentPosition);
        
        sum(sum((xHatTmp*rotationMatrix < 0) | (xHatTmp*rotationMatrix > 1)))
        
        % Rotation
        xHatTmp = xHatTmp*rotationMatrix;
        nuHat = nuHat*rotationMatrix;
        
        sum(sum((xHatTmp < 0) | (xHatTmp > 1)))
        sum(sum((xHatTmp*xHatTmp' < 0) | (xHatTmp*xHatTmp' > 1)))
        
        nv0 = (xHatTmp < 0);
        xHatTmp(nv0) = 0.0;
        
        nv1 = (sum(xHatTmp.^2, 2) > 1);
        xHatTmp(nv1, :) = xHatTmp(nv1, :)./...
            repmat(sum(xHatTmp(nv1, :).^2, 2), 1, dimLatentPosition);
        
        
        %% --- Pre-projection Method 2 ---
%         centerPoint = 1/(1 + sqrt(dimLatentPosition))*...
%             ones(nVertex, dimLatentPosition);
%         [~, xHatTmp1, transMatrix] = procrustes_revised(centerPoint, xHatTmp, ...
%             'scaling', false, 'reflection', false);
%         % norm(transMatrix.b*xHatTmp*transMatrix.T + transMatrix.c - xHatTmp1)
%         sum(sum((transMatrix.b*xHatTmp*transMatrix.T < 0) | (transMatrix.b*xHatTmp*transMatrix.T > 1)))
%         sum(sum((transMatrix.b*xHatTmp*transMatrix.T + transMatrix.c < 0) | ... 
%         (transMatrix.b*xHatTmp*transMatrix.T + transMatrix.c > 1)))
        
        fValueOld = objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
            adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
        
        while (~hasConverge) && (iter < maxIter)
            iter = iter + 1

            % Take 1 step of xHatTmp based on current nuHat, tauHat.
            [~, grad] = objectivefun(...
                reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, ...
                tauHat);
            
            grad = grad/max(max(abs(grad)));
            
            step = 1;
            xHatCurrent = xHatTmp - step*grad;
            pHatCurrent = xHatCurrent * xHatCurrent';
            while (any(any(xHatCurrent < zeros(nVertex, dimLatentPosition))) || ...
                    any(any(xHatCurrent > ones(nVertex, dimLatentPosition))) || ...
                    any(any(pHatCurrent < zeros(nVertex, nVertex))) || ...
                    any(any(pHatCurrent > ones(nVertex, nVertex))))
                step = step/2
                xHatCurrent = xHatTmp - step*grad;
                pHatCurrent = xHatCurrent * xHatCurrent';
            end
            xHatTmp = xHatTmp - step*grad;
            fValueNew = objectivefun(...
                reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, ...
                tauHat);
            
            % Convergence Checking
            if (abs(fValueNew - fValueOld) < tol*fValueOld) && (iter > 1)
                hasConverge = 1;
            end
            
            [tauHat, nuHat] = kmeans(xHatTmp, nBlock);
            tauHat = tauHat';
            
            fValueTmp = objectivefun(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                adjMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
            if (fValueTmp < fValueBest) || (iter == 1)
                fValueBest = fValueTmp;
                tauBest = tauHat;
            end
            
            fValueOld = fValueNew;
        end
        
        %% --- Save Results ---
        if (hasConverge == 1)
            errorRateOpti = errorratecalculator(tauStar, tauHat, ...
                nVertex, nBlock);
            errorRateOptiBest = errorratecalculator(tauStar, tauBest, ...
                nVertex, nBlock);
            parsave(saveFile, errorRateASGE, errorRateOpti, errorRateOptiBest);
        end
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

