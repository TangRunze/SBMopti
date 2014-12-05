close all
clear all

%% --- Minh's Model SBM ---

% Parameters setting

nBlock = 2;
dimLatentPosition = nBlock;
diag = 0.42;
offdiag = 0.42;
gStart = 1;
gEnd = 1000;

rho = [0.6 0.4];
B = [0.42 0.42; 0.42 0.5];

nVec = [100 250 500 750 1000];
nN = length(nVec);

rVec = [-1];
nR = length(rVec);

lambdaVec = [0.001 0.01 0.1 0.9 0.99 0.999];
nLambda = length(lambdaVec);

errorRateMinhSBMOptiVec = zeros(nR, nN);
errorRateMinhSBMMean = zeros(nR, nN);
errorRateMinhSBMCI = zeros(nR, nN, 2);
errorRateMinhSBMASGEVec = zeros(nR, nN);
errorRateMinhSBMCompare = zeros(nR, nN);
pValueMinhSBM = zeros(nR, nN);

for iN = 1:nN;
    nVertex = nVec(iN);
    
    % Go over different r.
    for iR = 1:nR
        r = rVec(iR);
        LOptiTmp = [];
        LASGETmp = [];
        % Go over all Monte Carlo replicates.
        for iGraph = gStart:gEnd
            LOptiVec = [];
            LASGEVec = [];
            loglikVec = [];
            % Go over different lambda.
            for iLambda = 1:nLambda
                lambda = lambdaVec(iLambda);
                if exist(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                        num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                        num2str(offdiag) '-r' num2str(r) '-lambda' ...
                        num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
                    load(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                        num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                        num2str(offdiag) '-r' num2str(r) '-lambda' ...
                        num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
                    LOptiVec = [LOptiVec, errorRateOptiBest];
                    LASGEVec = [LASGEVec, errorRateASGE];
                    loglikVec = [loglikVec, loglik];
                end
            end
            if (~isempty(LOptiVec))
                [~, ind] = max(loglikVec);
                LOptiTmp = [LOptiTmp LOptiVec(ind)];
                LASGETmp = [LASGETmp LASGEVec(ind)];
            end
        end
        errorRateMinhSBMOptiVec(iR, iN) = median(LOptiTmp);
        errorRateMinhSBMMean(iR, iN) = mean(LOptiTmp);
        errorRateMinhSBMCI(iR, iN, :) = [errorRateMinhSBMMean(iR, iN) - ...
            1.96*std(LOptiTmp)/sqrt(length(LOptiTmp)); ...
            errorRateMinhSBMMean(iR, iN) + ...
            1.96*std(LOptiTmp)/sqrt(length(LOptiTmp))];
        errorRateMinhSBMASGEVec(iR, iN) = median(LASGETmp);
        errorRateMinhSBMCompare(iR, iN) = sum(LOptiTmp < LASGETmp)/length(LOptiTmp);
        
        pValueMinhSBM(iR, iN) = 1 - binocdf(sum(LOptiTmp < LASGETmp) - 1, ...
            length(LOptiTmp) - sum(LOptiTmp == LASGETmp), 0.5);
    end
end

errorRateMinhSBMOptiVec
errorRateMinhSBMASGEVec
pValueMinhSBM

errorRateMinhSBMMean
errorRateMinhSBMCI

errorRateMinhSBMLambdaVec = zeros(nLambda, nN);
r = rVec(1);

for iN = 1:nN;
    nVertex = nVec(iN);
    
    % Go over different lambda.
    for iLambda = 1:nLambda
        lambda = lambdaVec(iLambda);
        LTmp = [];
        % Go over all Monte Carlo replicates.
        for iGraph = gStart:gEnd
            if exist(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                    num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                    num2str(offdiag) '-r' num2str(r) '-lambda' ...
                    num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
                load(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                    num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                    num2str(offdiag) '-r' num2str(r) '-lambda' ...
                    num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
                LTmp = [LTmp, errorRateOptiBest];
            end
        end
        errorRateMinhSBMLambdaVec(iLambda, iN) = mean(LTmp);
    end
end

errorRateMinhSBMLambdaVec

%% --- Minh's Model Dirichlet ---

nBlock = 2;
dimLatentPosition = nBlock;
diag = 0.42;
offdiag = 0.42;
gStart = 1;
gEnd = 1000;

rho = [0.6 0.4];
B = [0.42 0.42; 0.42 0.5];

nVec = [500];
nN = length(nVec);

rVec = [100];
nR = length(rVec);

lambdaVec = [0.001 0.01 0.1 0.9 0.99 0.999];
nLambda = length(lambdaVec);

errorRateMinhDirOptiVec = zeros(nR, nN);
errorRateMinhDirOptiMean = zeros(nR, nN);
errorRateMinhDirASGEVec = zeros(nR, nN);
errorRateMinhDirCompare = zeros(nR, nN);
pValueMinhDir = zeros(nR, nN);

for iN = 1:nN;
    nVertex = nVec(iN);
    
    % Go over different r.
    for iR = 1:nR
        r = rVec(iR);
        LOptiTmp = [];
        LASGETmp = [];
        % Go over all Monte Carlo replicates.
        for iGraph = gStart:gEnd
            LOptiVec = [];
            LASGEVec = [];
            loglikVec = [];
            % Go over different lambda.
            for iLambda = 1:nLambda
                lambda = lambdaVec(iLambda);
                if exist(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                        num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                        num2str(offdiag) '-r' num2str(r) '-lambda' ...
                        num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
                    load(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                        num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                        num2str(offdiag) '-r' num2str(r) '-lambda' ...
                        num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
                    LOptiVec = [LOptiVec, errorRateOptiBest];
                    LASGEVec = [LASGEVec, errorRateASGE];
                    loglikVec = [loglikVec, loglik];
                end
            end
            if (~isempty(LOptiVec))
                [~, ind] = max(loglikVec);
                LOptiTmp = [LOptiTmp LOptiVec(ind)];
                LASGETmp = [LASGETmp LASGEVec(ind)];
            end
        end
        errorRateMinhDirOptiVec(iR, iN) = median(LOptiTmp);
        errorRateMinhDirOptiMean(iR, iN) = mean(LOptiTmp);
        errorRateMinhDirASGEVec(iR, iN) = median(LASGETmp);
        errorRateMinhDirCompare(iR, iN) = sum(LOptiTmp < LASGETmp)/length(LOptiTmp);
        
        pValueMinhDir(iR, iN) = 1 - binocdf(sum(LOptiTmp < LASGETmp) - 1, ...
            length(LOptiTmp) - sum(LOptiTmp == LASGETmp), 0.5);
    end
end

errorRateMinhDirOptiVec
errorRateMinhDirASGEVec
pValueMinhDir

errorRateMinhDirOptiMean

%% --- Vince's Model SBM ---

nBlock = 3;
dimLatentPosition = nBlock;
muB = 0.5;
epsilonInB = 0.1;
diag = 0.6;
offdiag = 0.4;
gStart = 1;
gEnd = 1000;

rho = repmat(1/nBlock, 1, nBlock);
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

nVec = [150 300];
nN = length(nVec);

rVec = [-1];
nR = length(rVec);

lambdaVec = [0.001 0.01 0.1 0.9 0.99 0.999];
nLambda = length(lambdaVec);

errorRateVinceSBMOptiVec = zeros(nR, nN);
errorRateVinceSBMMean = zeros(nR, nN);
errorRateVinceSBMCI = zeros(nR, nN, 2);
errorRateVinceSBMASGEVec = zeros(nR, nN);
errorRateVinceSBMCompare = zeros(nR, nN);
pValueVinceSBM = zeros(nR, nN);

for iN = 1:nN;
    nVertex = nVec(iN);
    
    % Go over different r.
    for iR = 1:nR
        r = rVec(iR);
        LOptiTmp = [];
        LASGETmp = [];
        % Go over all Monte Carlo replicates.
        for iGraph = gStart:gEnd
            LOptiVec = [];
            LASGEVec = [];
            loglikVec = [];
            % Go over different lambda.
            for iLambda = 1:nLambda
                lambda = lambdaVec(iLambda);
                if exist(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                        num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                        num2str(offdiag) '-r' num2str(r) '-lambda' ...
                        num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
                    load(['./results/EB4SBMextended/results-SBMopti-Block-sim-dir-n' ...
                        num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                        num2str(offdiag) '-r' num2str(r) '-lambda' ...
                        num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
                    LOptiVec = [LOptiVec, errorRateOptiBest];
                    LASGEVec = [LASGEVec, errorRateASGE];
                    loglikVec = [loglikVec, loglik];
                end
            end
            if (~isempty(LOptiVec))
                [~, ind] = max(loglikVec);
                LOptiTmp = [LOptiTmp LOptiVec(ind)];
                LASGETmp = [LASGETmp LASGEVec(ind)];
            end
        end
        errorRateVinceSBMOptiVec(iR, iN) = median(LOptiTmp);
        errorRateVinceSBMMean(iR, iN) = mean(LOptiTmp);
        errorRateVinceSBMCI(iR, iN, :) = [errorRateVinceSBMMean(iR, iN) - ...
            1.96*std(LOptiTmp)/sqrt(length(LOptiTmp)); ...
            errorRateVinceSBMMean(iR, iN) + ...
            1.96*std(LOptiTmp)/sqrt(length(LOptiTmp))];
        errorRateVinceSBMASGEVec(iR, iN) = median(LASGETmp);
        errorRateVinceSBMCompare(iR, iN) = sum(LOptiTmp < LASGETmp)/length(LOptiTmp);
        
        pValueVinceSBM(iR, iN) = 1 - binocdf(sum(LOptiTmp < LASGETmp) - 1, ...
            length(LOptiTmp) - sum(LOptiTmp == LASGETmp), 0.5);
    end
end

errorRateVinceSBMOptiVec
errorRateVinceSBMASGEVec
pValueVinceSBM

errorRateVinceSBMMean
errorRateVinceSBMCI