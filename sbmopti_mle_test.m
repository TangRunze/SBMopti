

function sbmopti_mle_test(r)

%% --- Parameter Setting ---

nVertex = 150;
nBlock = 3;
dimLatentPosition = nBlock;
rho = repmat(1/nBlock, 1, nBlock);
maxIter = 100;
nCore = 12;
muB = 0.5;
epsilonInB = 0.1;
% r = -1;
tol = 1e-4;

gStart = 1;
gEnd = 1000;

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
options = optimoptions('fmincon', 'TolX', 1e-6, ...
    'MaxIter', 5000, 'MaxFunEvals', 5000, ...
    'Algorithm', 'interior-point'); %, 'GradObj', 'on');
% Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);

%% --- Parallel Computing ---
delete(gcp('nocreate'))
parpool(nCore);

lambdaVec = [0 0.01 0.1];

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
            
            if (lambda == 0)
                hasConverge = 1;
                tauHat = tauHat0;
                tauBest = tauHat;
                xHatTmp = xHat0;
                xBest = xHatTmp;
                nuHat = nuHat0;
                nuBest = nuHat;
                fValueBest = objectivefun_std(reshape(xHatTmp, 1, ...
                    nVertex*dimLatentPosition), adjMatrix, nuHat, ...
                    lambda, nVertex, dimLatentPosition, tauHat);
                fValueNew = fValueBest;

            else
                
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

% %% --- Close Parallel Computing ---
delete(gcp('nocreate'))




%% Analysis

% ind = [];
% for iGraph = gStart:gEnd
%     flag = true;
%     for lambda = lambdaVec
%         if ~exist(['./results/results-SBMopti-Block-sim-dir-n' ...
%                 num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
%                 num2str(B(1, 2)) '-r' num2str(r) '-lambda' num2str(lambda) ...
%                 '-adjMatrix' num2str(iGraph) '.mat'])
%             flag = false;
% 
%         end
%     end
%     if (flag == true)
%         ind = [ind iGraph];
%     end
% end
% 
% maxIter = length(ind);
% loglikVec = zeros(maxIter, length(lambdaVec));
% 
% for iInd = 1:maxIter
%     for iLambda = 1:length(lambdaVec)
%         lambda = lambdaVec(iLambda);
%         load(['./results/results-SBMopti-Block-sim-dir-n' ...
%             num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
%             num2str(B(1, 2)) '-r' num2str(r) '-lambda' num2str(lambda) ...
%             '-adjMatrix' num2str(ind(iInd)) '.mat']);
%         loglikVec(iInd, iLambda) = loglik;
%         errorRateVec(iInd, iLambda) = errorRateOptiBest;
%     end
% end
% 
% sum(loglikVec(:, 1) < loglikVec(:, 2))/maxIter
% sum(errorRateVec(:, 1) > errorRateVec(:, 2))/maxIter

