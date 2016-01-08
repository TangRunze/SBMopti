close all
clear all

%% --- Simulation ---

% pathExtra = '/';
pathExtra = '/sim3/';


% Parameter Setting
nVertex = 150;
gStart = 1;
gEnd = 200;

B = [0.4, 0.2, 0.2; 0.2, 0.4, 0.2; 0.2, 0.2, 0.4];

lambdaVec = 0.999;

dimLatentPosition = 3;
isGMM = 0;

r = -1;

ind = [];
for iGraph = gStart:gEnd
    iLambda = 1;
    lambda = lambdaVec(iLambda);
    while (iLambda <= length(lambdaVec)) && ...
            (exist(['./results' pathExtra 'results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-GMM' num2str(isGMM) '-r' num2str(r) ...
            '-lambda' num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']) == 0)
        iLambda = iLambda + 1;
        if iLambda <= length(lambdaVec)
            lambda = lambdaVec(iLambda);
        end
    end
    if (iLambda <= length(lambdaVec))
        ind = [ind iGraph];
    end
end

maxIter = length(ind);

errorASGEVec = zeros(1, maxIter);
errorOptiVec = zeros(1, maxIter);
loglikASGEVec = zeros(1, maxIter);
loglikOptiVec = zeros(1, maxIter);

for iInd = 1:maxIter
    iGraph = ind(iInd);
    if exist(['./results' pathExtra 'results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-GMM' num2str(isGMM) '-r' num2str(r) ...
            '-lambda' num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
        load(['./results' pathExtra 'results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
            num2str(B(1, 2)) '-GMM' num2str(isGMM) '-r' num2str(r) ...
            '-lambda' num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
        errorASGEVec(iInd) = errorRateASGE;
        errorOptiVec(iInd) = errorRateOpti;
        loglikASGEVec(iInd) = loglikASGE;
        loglikOptiVec(iInd) = loglikOpti;
    end
end

errorMeanASGE = mean(errorASGEVec);
errorCIASGE = [errorMeanASGE - 1.96*std(errorASGEVec)/sqrt(maxIter), ...
    errorMeanASGE + 1.96*std(errorASGEVec)/sqrt(maxIter)];
errorMedianASGE = median(errorASGEVec);

errorMeanOpti = mean(errorOptiVec);
errorCIOpti = [errorMeanOpti - 1.96*std(errorOptiVec)/sqrt(maxIter), ...
    errorMeanOpti + 1.96*std(errorOptiVec)/sqrt(maxIter)];
errorMedianOpti = median(errorOptiVec);

loglikMeanASGE = mean(loglikASGEVec);
loglikMeanOpti = mean(loglikOptiVec);

[errorMeanASGE, errorMeanOpti]
% [loglikMeanASGE, loglikMeanOpti]

%% --- Paired Test ---

% % Alternative hypothesis: Opti < ASGE.
% 
% errorPair = zeros(2, maxIter);
% errorPair(1, :) = errorOptiVec;
% errorPair(2, :) = errorASGEVec;
% 
% % 1-sided sign-test
% tmpStats = sum(errorPair(1, :) < errorPair(2, :));
% pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
%     sum(errorPair(1, :) == errorPair(2, :)), 0.5)



% errorPair = zeros(2, maxIter);
% errorPair(1, :) = loglikOptiVec;
% errorPair(2, :) = loglikASGEVec;
% 
% % 1-sided sign-test
% tmpStats = sum(errorPair(1, :) > errorPair(2, :));
% pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
%     sum(errorPair(1, :) == errorPair(2, :)), 0.5)
% 
% 
% 
% errorPair = zeros(2, maxIter);
% errorPair(1, :) = loglikOpticVec;
% errorPair(2, :) = loglikASGEcVec;
% 
% % 1-sided sign-test
% tmpStats = sum(errorPair(1, :) > errorPair(2, :));
% pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
%     sum(errorPair(1, :) == errorPair(2, :)), 0.5)

% % Plot
% figure;
% tmpStats = round((errorPair(1, :) - errorPair(2, :))*nVertex);
% hist(tmpStats, -80:10:80)
% 
% xlabel('Optimization Error - ASGE Error');
% ylabel('Frequency')





