close all
clear all

%% --- Simulation ---

% Parameter Setting
nVertex = 150;
epsilon = 0.1;
gStart = 1;
gEnd = 1000;

% lambdaRange = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 10, 20, 50, 100];
lambdaRange = (1:10)/10;

% ind = [];
% for iGraph = gStart:gEnd
%     if exist(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
%             num2str(epsilon) '-lambda' num2str(lambda) ...
%             '-graph' num2str(iGraph) '.mat'])
%         ind = [ind iGraph];
%     end
% end

ind = gStart:gEnd;

maxIter = length(ind);

errorASGE = zeros(length(lambdaRange), maxIter);
errorOpti = zeros(length(lambdaRange), maxIter);
errorOptiBest = zeros(length(lambdaRange), maxIter);
fValue = zeros(length(lambdaRange), maxIter);
fValueBest = zeros(length(lambdaRange), maxIter);

errorMedianASGE = zeros(1, length(lambdaRange));
errorMedianOpti = zeros(1, length(lambdaRange));
errorMedianOptiBest = zeros(1, length(lambdaRange));
medianfValue = zeros(1, length(lambdaRange));
medianfValueBest = zeros(1, length(lambdaRange));

for iLambda = 1:length(lambdaRange)
    lambda = lambdaRange(iLambda);
    for iInd = 1:maxIter
        load(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-lambda' num2str(lambda) ...
            '-graph' num2str(ind(iInd)) '.mat']);
        errorASGE(iLambda, iInd) = errorRateASGE;
        errorOpti(iLambda, iInd) = errorRateOpti;
        errorOptiBest(iLambda, iInd) = errorRateOptiBest;
        fValue(iLambda, iInd) = fValueNew;
        fValueBest(iLambda, iInd) = fValueBest;
    end
    
    % errorMeanASGE = mean(errorASGE);
    % errorCIASGE = [errorMeanASGE - 1.96*std(errorASGE)/sqrt(maxIter), ...
    %     errorMeanASGE + 1.96*std(errorASGE)/sqrt(maxIter)];
    errorMedianASGE(iLambda) = median(errorASGE(iLambda, :));
    
    % errorMeanOpti = mean(errorOpti);
    % errorCIOpti = [errorMeanOpti - 1.96*std(errorOpti)/sqrt(maxIter), ...
    %     errorMeanOpti + 1.96*std(errorOpti)/sqrt(maxIter)];
    errorMedianOpti(iLambda) = median(errorOpti(iLambda, :));
    
    errorMedianOptiBest(iLambda) = median(errorOptiBest(iLambda, :));
    
    medianfValue(iLambda) = median(fValue(iLambda, :));
    
    medianfValueBest(iLambda) = median(fValueBest(iLambda, :));
    
end

%% --- Plot ---

figure;
plot(lambdaRange, medianfValue, 'go');
hold on;
plot(lambdaRange, medianfValue, 'b');
xlabel('lambda');
ylabel('objective function value');
hold off;

figure;
plot(lambdaRange, errorMedianOpti, 'go');
hold on;
plot(lambdaRange, errorMedianOpti, 'r');
xlabel('lambda');
ylabel('Lhat');
hold off;



% figure;
% plot(lambdaRange(1:8), medianfValue(1:8), 'go');
% hold on;
% plot(lambdaRange(1:8), medianfValue(1:8), 'b');
% xlabel('lambda');
% ylabel('objective function value');
% hold off;
% 
% figure;
% plot(lambdaRange(1:8), errorMedianOpti(1:8), 'go');
% hold on;
% plot(lambdaRange(1:8), errorMedianOpti(1:8), 'r');
% xlabel('lambda');
% ylabel('Lhat');
% hold off;

%% --- Paired Test ---

% Alternative hypothesis: pair1 < pair2.

lambda1 = 0.1;
lambda2 = 0.2;

indPair = [];
for iGraph = gStart:gEnd
    if exist(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-lambda' num2str(lambda1) ...
            '-graph' num2str(ind(iGraph)) '.mat']) && ...
            exist(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-lambda' num2str(lambda2) ...
            '-graph' num2str(ind(iGraph)) '.mat'])
        indPair = [indPair iGraph];
    end
end

maxIterPair = length(indPair);

fValue = zeros(2, maxIterPair);

for iIndPair = 1:maxIterPair
    load(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-lambda' num2str(lambda1) ...
            '-graph' num2str(ind(iIndPair)) '.mat']);
    fValue(1, iIndPair) = fValueNew;
    load(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-lambda' num2str(lambda2) ...
            '-graph' num2str(ind(iIndPair)) '.mat']);
    fValue(2, iIndPair) = fValueNew;
end


% 1-sided sign-test
tmpStats = sum(fValue(1, :) < fValue(2, :));
pValue = 1 - binocdf(tmpStats - 1, maxIterPair, 0.5)







tmpStats = round((fValue(1, :) - fValue(2, :))*nVertex);
hist(tmpStats, -60:1:60)



[n, xout] = hist(tmpStats);
bar(xout, n, 1)




