close all
clear all

%% --- Simulation ---

% Parameter Setting
nVertex = 300;
gStart = 1;
gEnd = 1000;

ind = [];
for iGraph = gStart:gEnd
    if exist(['./results/results-SBMopti-real-graph' num2str(iGraph) '.mat'])
        ind = [ind iGraph];
    end
end

maxIter = length(ind);

errorASGE = zeros(1, maxIter);
errorOpti = zeros(1, maxIter);

for iInd = 1:maxIter
    load(['./results/results-SBMopti-real-graph' num2str(ind(iInd)) '.mat']);
    errorASGE(iInd) = errorRateASGE;
    errorOpti(iInd) = errorRateOpti;
end

errorMeanASGE = mean(errorASGE);
errorCIASGE = [errorMeanASGE - 1.96*std(errorASGE)/sqrt(maxIter), ...
    errorMeanASGE + 1.96*std(errorASGE)/sqrt(maxIter)];
errorMedianASGE = median(errorASGE);

errorMeanOpti = mean(errorOpti);
errorCIOpti = [errorMeanOpti - 1.96*std(errorOpti)/sqrt(maxIter), ...
    errorMeanOpti + 1.96*std(errorOpti)/sqrt(maxIter)];
errorMedianOpti = median(errorOpti);


%% --- Plot ---

binStart = 0;
binEnd = 0.6;
nBin = 20;
myBins = linspace(binStart + (binEnd - binStart)/(nBin - 1), ...
    binEnd - (binEnd - binStart)/(nBin - 1), nBin);

hist(errorOpti, myBins);
hold on;
hist(errorASGE, myBins);

h = findobj(gca, 'Type', 'patch');

set(h(1), 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75);
set(h(2), 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75);

legend('Optimization', 'ASGE');
legend boxoff

xlabel('Error Rate');
ylabel('Frequency')

set(gca,'box','off');

hold off

%% --- Paired Test ---

% Alternative hypothesis: Opti < ASGE.

errorPair = zeros(2, maxIter);
errorPair(1, :) = errorOpti;
errorPair(2, :) = errorASGE;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) < errorPair(2, :));
pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) < errorPair(2, :)), 0.5)

% Plot
figure;
tmpStats = round((errorPair(1, :) - errorPair(2, :))*nVertex);
hist(tmpStats, -80:10:80)

xlabel('Optimization Error - ASGE Error');
ylabel('Frequency')





