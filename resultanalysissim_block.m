close all
clear all

%% --- Simulation ---

% Parameter Setting
nVertex = 150;
epsilon = 0.02;
gStart = 1;
gEnd = 1000;

delta = 0.5;

% lambdaVec = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 10, 20, 50, 100];
% lambdaVec = [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000];
lambdaVec = [0.001 0.01 0.1 0.5 0.9 0.99 0.999];

errorASGE = [];
errorOpti = [];
errorOptiBest = [];
fValue = [];
fValueBest = [];

errorMedianASGE = zeros(1, length(lambdaVec));
errorMedianOpti = zeros(1, length(lambdaVec));
errorMedianOptiBest = zeros(1, length(lambdaVec));
medianfValue = zeros(1, length(lambdaVec));
medianfValueBest = zeros(1, length(lambdaVec));

for iLambda = 1:length(lambdaVec)
    lambda = lambdaVec(iLambda);
    for iGraph = gStart:gEnd
        if exist(['./results/results-SBMopti-Block-sim-n' ...
                num2str(nVertex) '-eps' num2str(epsilon) '-delta' ...
                num2str(delta) '-lambda' num2str(lambda) '-pmatrix' ...
                num2str(iGraph) '.mat'])
            load(['./results/results-SBMopti-Block-sim-n' ...
                num2str(nVertex) '-eps' num2str(epsilon) '-delta' ...
                num2str(delta) '-lambda' num2str(lambda) '-pmatrix' ...
                num2str(iGraph) '.mat']);
            errorASGE = [errorASGE errorRateASGE];
            % errorOpti = [errorOpti errorRateOpti];
            errorOptiBest = [errorOptiBest errorRateOptiBest];
            % fValue = [fValue fValueNew];
            fValueBest = [fValueBest fValueBest];
        end
    end
    
    % errorMeanASGE = mean(errorASGE);
    % errorCIASGE = [errorMeanASGE - 1.96*std(errorASGE)/sqrt(maxIter), ...
    %     errorMeanASGE + 1.96*std(errorASGE)/sqrt(maxIter)];
    errorMedianASGE(iLambda) = median(errorASGE);

    % errorMeanOpti = mean(errorOpti);
    % errorCIOpti = [errorMeanOpti - 1.96*std(errorOpti)/sqrt(maxIter), ...
    %     errorMeanOpti + 1.96*std(errorOpti)/sqrt(maxIter)];
    errorMedianOpti(iLambda) = median(errorOpti);

    errorMedianOptiBest(iLambda) = median(errorOptiBest);

    medianfValue(iLambda) = median(fValue);

    medianfValueBest(iLambda) = median(fValueBest);

end

errorMedianASGE

errorMedianOptiBest

medianfValueBest


%% --- Blockness ---

% nVertex = 150;
% epsilon = 0.02;
% gStart = 1;
% gEnd = 1000;
% 
% f0Median = zeros(1, 9);
% f1Median = zeros(1, 9);
% fIntersectMedian = zeros(1, 9);
% f0Mean = zeros(1, 9);
% f1Mean = zeros(1, 9);
% fIntersectMean = zeros(1, 9);
% f0CI = zeros(2, 9);
% f1CI = zeros(2, 9);
% fIntersectCI = zeros(2, 9);
% 
% for iDelta = 1:9
%     delta = iDelta/10;
%     
%     f0Tmp = [];
%     f1Tmp = [];
%     fIntersectTmp = [];
%     for iGraph = gStart:gEnd
%         if exist(['./results/results-SBMopti-Block-sim-n' ...
%                 num2str(nVertex) '-eps' num2str(epsilon) '-delta' ...
%                 num2str(delta) '-lambda' num2str(0) '-pmatrix' ...
%                 num2str(iGraph) '.mat']) && ...
%                 exist(['./results/results-SBMopti-Block-sim-n' ...
%                 num2str(nVertex) '-eps' num2str(epsilon) '-delta' ...
%                 num2str(delta) '-lambda' num2str(1) '-pmatrix' ...
%                 num2str(iGraph) '.mat'])
%             load(['./results/results-SBMopti-Block-sim-n' ...
%                 num2str(nVertex) '-eps' num2str(epsilon) '-delta' ...
%                 num2str(delta) '-lambda' num2str(0) '-pmatrix' ...
%                 num2str(iGraph) '.mat']);
%             f1Tmp = [f1Tmp fValue1];
%             load(['./results/results-SBMopti-Block-sim-n' ...
%                 num2str(nVertex) '-eps' num2str(epsilon) '-delta' ...
%                 num2str(delta) '-lambda' num2str(1) '-pmatrix' ...
%                 num2str(iGraph) '.mat']);
%             f0Tmp = [f0Tmp fValue0];
%             fIntersectTmp = [fIntersectTmp f0Tmp(end)*f1Tmp(end)/(f0Tmp(end)+f1Tmp(end))];
%         end
%     end
%     
%     f0Median(iDelta) = median(f0Tmp);
%     f1Median(iDelta) = median(f1Tmp);
%     fIntersectMedian(iDelta) = median(fIntersectTmp);
%     
%     f0Mean(iDelta) = mean(f0Tmp);
%     f1Mean(iDelta) = mean(f1Tmp);
%     fIntersectMean(iDelta) = mean(fIntersectTmp);
%     
%     f0CI(:, iDelta) = [f0Mean(iDelta) - 1.96*std(f0Tmp)/sqrt(length(f0Tmp)); ...
%         f0Mean(iDelta) + 1.96*std(f0Tmp)/sqrt(length(f0Tmp))];
%     f1CI(:, iDelta) = [f1Mean(iDelta) - 1.96*std(f1Tmp)/sqrt(length(f1Tmp)); ...
%         f1Mean(iDelta) + 1.96*std(f1Tmp)/sqrt(length(f1Tmp))];
%     fIntersectCI(:, iDelta) = [fIntersectMean(iDelta) - 1.96*std(fIntersectTmp)/sqrt(length(fIntersectTmp)); ...
%         fIntersectMean(iDelta) + 1.96*std(fIntersectTmp)/sqrt(length(fIntersectTmp))];
%     
% end
% 
% f0Median
% f1Median
% fIntersectMedian
% 
% f0CI
% f1CI
% fIntersectCI


%% --- Plot ---

% figure;
% plot(lambdaVec, medianfValue, 'go');
% hold on;
% plot(lambdaVec, medianfValue, 'b');
% xlabel('lambda');
% ylabel('objective function value');
% hold off;
% 
% figure;
% plot(lambdaVec, errorMedianOpti, 'go');
% hold on;
% plot(lambdaVec, errorMedianOpti, 'r');
% xlabel('lambda');
% ylabel('Lhat');
% hold off;



%% --- Paired Test ---
% 
% % Alternative hypothesis: pair1 < pair2.
% 
% lambda1 = 0.1;
% lambda2 = 0.2;
% 
% indPair = [];
% for iGraph = gStart:gEnd
%     if exist(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
%             num2str(epsilon) '-lambda' num2str(lambda1) ...
%             '-graph' num2str(ind(iGraph)) '.mat']) && ...
%             exist(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
%             num2str(epsilon) '-lambda' num2str(lambda2) ...
%             '-graph' num2str(ind(iGraph)) '.mat'])
%         indPair = [indPair iGraph];
%     end
% end
% 
% maxIterPair = length(indPair);
% 
% fValue = zeros(2, maxIterPair);
% 
% for iIndPair = 1:maxIterPair
%     load(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
%         num2str(epsilon) '-lambda' num2str(lambda1) ...
%         '-graph' num2str(ind(iIndPair)) '.mat']);
%     fValue(1, iIndPair) = fValueNew;
%     load(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
%         num2str(epsilon) '-lambda' num2str(lambda2) ...
%         '-graph' num2str(ind(iIndPair)) '.mat']);
%     fValue(2, iIndPair) = fValueNew;
% end
% 
% 
% % 1-sided sign-test
% tmpStats = sum(fValue(1, :) < fValue(2, :));
% pValue = 1 - binocdf(tmpStats - 1, maxIterPair, 0.5)
% 
% 
% 
% 
% 
% 
% 
% tmpStats = round((fValue(1, :) - fValue(2, :))*nVertex);
% hist(tmpStats, -60:1:60)
% 
% 
% 
% [n, xout] = hist(tmpStats);
% bar(xout, n, 1)



