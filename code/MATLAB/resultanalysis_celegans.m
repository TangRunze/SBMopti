close all
clear all

%% --- Simulation ---

% Parameter Setting
nVertex = 300;
gStart = 1;
gEnd = 200;

lambdaVec = 0.5;
lambda = 0.5;

dimLatentPosition = 6;
isGMM = 1;

ind = [];
for iGraph = gStart:gEnd
    iLambda = 1;
    lambda = lambdaVec(iLambda);
    while (iLambda <= length(lambdaVec)) && ...
            (exist(['./results/results-SBMopti-celegans-graph' ...
            num2str(iGraph) '-dim' num2str(dimLatentPosition) ...
            '-GMM' num2str(isGMM) '-lambda' num2str(lambda) '.mat']) == 0)
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
loglikASGEcVec = zeros(1, maxIter);
loglikOpticVec = zeros(1, maxIter);

nv = [];

for iInd = 1:maxIter
    iGraph = ind(iInd);
    if exist(['./results/results-SBMopti-celegans-graph' ...
            num2str(iGraph) '-dim' num2str(dimLatentPosition) ...
            '-GMM' num2str(isGMM) '-lambda' num2str(lambda) '.mat'])
        load(['./results/results-SBMopti-celegans-graph' ...
            num2str(iGraph) '-dim' num2str(dimLatentPosition) ...
            '-GMM' num2str(isGMM) '-lambda' num2str(lambda) '.mat']);
        errorASGEVec(iInd) = errorRateASGE;
        errorOptiVec(iInd) = errorRateOpti;
        loglikASGEVec(iInd) = loglikASGE;
        loglikOptiVec(iInd) = loglikOpti;
        loglikASGEcVec(iInd) = loglikASGEc;
        loglikOpticVec(iInd) = loglikOptic;
        if (isnan(loglikASGE) || isnan(loglikOpti) || ...
                isnan(loglikASGEc) || isnan(loglikOptic))
            nv = [nv, iInd];
        end
    end
end

ind(nv) = [];
errorASGEVec(nv) = [];
errorOptiVec(nv) = [];
loglikASGEVec(nv) = [];
loglikOptiVec(nv) = [];
loglikASGEcVec(nv) = [];
loglikOpticVec(nv) = [];
maxIter = length(ind);


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

loglikMeanASGEc = mean(loglikASGEcVec);
loglikMeanOptic = mean(loglikOpticVec);

%% --- Paired Test ---

% Alternative hypothesis: Opti < ASGE.

errorPair = zeros(2, maxIter);
errorPair(1, :) = errorOptiVec;
errorPair(2, :) = errorASGEVec;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) < errorPair(2, :));
pValue1 = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) == errorPair(2, :)), 0.5);






errorPair = zeros(2, maxIter);
errorPair(1, :) = loglikOpticVec;
errorPair(2, :) = loglikASGEcVec;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) > errorPair(2, :));
pValue2 = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) == errorPair(2, :)), 0.5);




errorPair = zeros(2, maxIter);
errorPair(1, :) = loglikOptiVec;
errorPair(2, :) = loglikASGEVec;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) > errorPair(2, :));
pValue3 = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) == errorPair(2, :)), 0.5);

% % Plot
% figure;
% tmpStats = round((errorPair(1, :) - errorPair(2, :))*nVertex);
% hist(tmpStats, -80:10:80)
% 
% xlabel('Optimization Error - ASGE Error');
% ylabel('Frequency')



%% Display result

[errorMeanASGE, errorMeanOpti]
pValue1
[loglikMeanASGEc, loglikMeanOptic]
pValue2
[loglikMeanASGE, loglikMeanOpti]
pValue3



