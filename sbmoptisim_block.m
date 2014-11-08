function [] = sbmoptisim_block(nVertex, nBlock, muB, epsilonInB, r, ...
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

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

%% --- Optimization Setting ---
options = optimoptions('fmincon', 'TolX', 1e-6, ...
    'MaxIter', 10000, 'MaxFunEvals', 10000, ...
    'Algorithm', 'interior-point'); %, 'GradObj', 'on');
% Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);

%% --- Parallel Computing ---
delete(gcp('nocreate'))
parpool(nCore);

lambdaVec = [0.001 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.999];

parfor iProbMatrix = gStart:gEnd
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [pMatrix, nuHat0, ~, tauHat0, ~] = ...
        datagenerator_block(nVertex, nBlock, dimLatentPosition, B, rho, ...
        tauStar, epsilonInB, r, iProbMatrix, projectoptions);
    xHat = asge(pMatrix, dimLatentPosition);
    
    % Pre-projection
    rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
        dimLatentPosition, xHat), ...
        reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
        [], [], [], [], - ones(dimLatentPosition^2, 1), ...
        ones(dimLatentPosition^2, 1), @(x) ...
        projectconditionfun(x, nVertex, dimLatentPosition, xHat), ...
        projectoptions);
    rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
        dimLatentPosition);
    
    % Rotate the latent positions
    xHat = xHat*rotationMatrix;
    nuHat0 = nuHat0*rotationMatrix;
    
    %% --- ASGE ---
    errorRateASGE = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);
    
    for lambda = lambdaVec
        
        saveFile = ['./results/results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-r' num2str(r) '-lambda' num2str(lambda) ...
            '-pmatrix' num2str(iProbMatrix) '.mat'];
        
        if exist(saveFile, 'file') == 0
            %% --- Solve Optimization Problem ---
            hasConverge = 0;
            iter = 0;
            fValueBest = 0;
            xHatTmp = xHat;
            nuHat = nuHat0;
            tauHat = tauHat0;
            
            fValueOld = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
            
            while (~hasConverge) && (iter < maxIter)
                iter = iter + 1
                
                % Find the best xHatTmp based on current nuHat, tauHat.
                xHatTmp = fmincon(@(x) objectivefun_std(x, ...
                    pMatrix, nuHat, lambda, nVertex, dimLatentPosition, ...
                    tauHat), reshape(xHatTmp, 1, nVertex*dimLatentPosition),...
                    [], [], [], [], zeros(dimLatentPosition*nVertex, 1), ...
                    ones(dimLatentPosition*nVertex, 1), @(x) ...
                    conditionfun(x, nVertex, dimLatentPosition), options);
                xHatTmp = reshape(xHatTmp, nVertex, dimLatentPosition);
                
                % Find the best nuHat, tauHat based on current xHatTmp.
                [tauHat, nuHat] = kmeans(xHatTmp, nBlock);
                tauHat = tauHat';
                
                fValueNew = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                    pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
                if (fValueNew < fValueBest) || (iter == 1)
                    fValueBest = fValueNew;
                    tauBest = tauHat;
                    xBest = xHatTmp;
                    nuBest = nuHat;
                    
                    fValue0 = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                        pMatrix, nuHat, 0, nVertex, dimLatentPosition, tauHat);
                    fValue1 = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
                        pMatrix, nuHat, 1, nVertex, dimLatentPosition, tauHat);
                    
                end
                
                % Convergence Checking
                if (abs(fValueNew - fValueOld) < tol*fValueOld) && (iter > 1)
                    hasConverge = 1;
                end
                
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
                    errorRateOptiBest, fValueNew, fValueBest, ...
                    xBest, nuBest, tauBest, fValue0, fValue1);
            end
        end
    end
    
    %% Consider lambda = 0
    lambda = 0;
    
    saveFile = ['./results/results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-r' num2str(r) '-lambda' num2str(lambda) ...
            '-pmatrix' num2str(iProbMatrix) '.mat']
    
    if exist(saveFile, 'file') == 0
        
        %% --- Solve Optimization Problem ---
        xHatTmp = xHat;
        
        [tauHat, nuHat] = kmeans(xHatTmp, nBlock);
        tauHat = tauHat';
        
        fValueNew = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
            pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
        fValueBest = fValueNew;
        
        fValue0 = fValueBest;
        fValue1 = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
            pMatrix, nuHat, 1, nVertex, dimLatentPosition, tauHat);
        
        %% --- Save Results ---
        errorRateOpti = errorratecalculator(tauStar, tauHat, ...
            nVertex, nBlock);
        errorRateOptiBest = errorRateOpti;
        parsave(saveFile, errorRateASGE, errorRateOpti, ...
            errorRateOptiBest, fValueNew, fValueBest, ...
            xHatTmp, nuHat, tauHat, fValue0, fValue1);
    end
    
    %% Consider lambda = 1
    lambda = 1;
    
    saveFile = ['./results/results-SBMopti-Block-sim-dir-n' ...
        num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
        num2str(B(1, 2)) '-r' num2str(r) '-lambda' num2str(lambda) ...
        '-pmatrix' num2str(iProbMatrix) '.mat'];
    
    if exist(saveFile, 'file') == 0
        
        %% --- Solve Optimization Problem ---
        nuHat = nuHat0;
        tauHat = tauHat0;
        
        % Find the best xHatTmp based on current nuHat, tauHat.
        nuHat = fmincon(@(x) objectivefun_std_lambda1(...
            pMatrix, x, dimLatentPosition, nBlock, tauHat), ...
            reshape(nuHat, 1, nBlock*dimLatentPosition),...
            [], [], [], [], zeros(dimLatentPosition*nBlock, 1), ...
            ones(dimLatentPosition*nBlock, 1), @(x) ...
            conditionfun(x, nBlock, dimLatentPosition), options);
        nuHat = reshape(nuHat, nBlock, dimLatentPosition);

        xHatTmp = nuHat(tauHat, :);
        
        fValueNew = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
            pMatrix, nuHat, lambda, nVertex, dimLatentPosition, tauHat);
        fValueBest = fValueNew;
        
        fValue1 = fValueBest;
        fValue0 = objectivefun_std(reshape(xHatTmp, 1, nVertex*dimLatentPosition), ...
            pMatrix, nuHat, 0, nVertex, dimLatentPosition, tauHat);
        
        %% --- Save Results ---
        errorRateOpti = errorratecalculator(tauStar, tauHat, ...
            nVertex, nBlock);
        errorRateOptiBest = errorRateOpti;
        parsave(saveFile, errorRateASGE, errorRateOpti, ...
            errorRateOptiBest, fValueNew, fValueBest, ...
            xHatTmp, nuHat, tauHat, fValue0, fValue1);
    end
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

