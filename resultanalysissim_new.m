close all
clear all

%% --- Simulation ---

% Parameter Setting
nVertex = 150;
epsilon = 0.1;
gStart = 1;
gEnd = 1000;

ind = [];
for iGraph = gStart:gEnd
    if exist(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-graph' num2str(iGraph) '.mat'])
        ind = [ind iGraph];
    end
end

maxIter = length(ind);

errorASGE = zeros(1, maxIter);
errorOpti = zeros(1, maxIter);
errorOptiBest = zeros(1, maxIter);

scatter1 = [];
scatter2 = [];

for iInd = 1:maxIter
    load(['./results/results-SBMopti-sim-n' num2str(nVertex) '-eps' ...
            num2str(epsilon) '-graph' num2str(ind(iInd)) '.mat']);
    errorASGE(iInd) = errorRateASGE;
    errorOpti(iInd) = errorRateOpti;
    errorOptiBest(iInd) = errorRateOptiBest;
    scatterdata = scatterdata(:, 1:(end - 1));
    scatter1 = [scatter1, scatterdata([1, 3], :)];
    scatter2 = [scatter2, scatterdata([2, 3], :)];
end

plot(scatter1(2, :), scatter1(1, :), 'o');
plot(scatter2(2, :), scatter1(1, :), 'o');

errorMeanASGE = mean(errorASGE);
errorCIASGE = [errorMeanASGE - 1.96*std(errorASGE)/sqrt(maxIter), ...
    errorMeanASGE + 1.96*std(errorASGE)/sqrt(maxIter)];
errorMedianASGE = median(errorASGE);

errorMeanOpti = mean(errorOpti);
errorCIOpti = [errorMeanOpti - 1.96*std(errorOpti)/sqrt(maxIter), ...
    errorMeanOpti + 1.96*std(errorOpti)/sqrt(maxIter)];
errorMedianOpti = median(errorOpti);

errorMeanOptiBest = mean(errorOptiBest);
errorCIOptiBest = [errorMeanOptiBest - 1.96*std(errorOptiBest)/sqrt(maxIter), ...
    errorMeanOptiBest + 1.96*std(errorOptiBest)/sqrt(maxIter)];
errorMedianOptiBest = median(errorOptiBest);

%% --- Plot ---

binStart = 0;
binEnd = 0.6;
nBin = 20;
myBins = linspace(binStart + (binEnd - binStart)/(nBin - 1), ...
    binEnd - (binEnd - binStart)/(nBin - 1), nBin);

hist(errorASGE, myBins);
hold on;
hist(errorOpti, myBins);
hold on;
hist(errorOptiBest, myBins);

h = findobj(gca, 'Type', 'patch');

set(h(1), 'FaceColor', 'g', 'EdgeColor', 'w', 'facealpha', 0.7);
set(h(2), 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.7);
set(h(3), 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.7);

legend('ASGE', 'Optimization', 'Optimization Best');
legend boxoff

xlabel('Error Rate');
ylabel('Frequency')

set(gca,'box','off');

hold off

%% --- Paired Test ---

errorPair = zeros(2, maxIter);

% Alternative hypothesis: Opti < ASGE.
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


% Alternative hypothesis: OptiBest < Opti.
errorPair(1, :) = errorOptiBest;
errorPair(2, :) = errorOpti;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) < errorPair(2, :));
pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) < errorPair(2, :)), 0.5)

% Plot
figure;
tmpStats = round((errorPair(1, :) - errorPair(2, :))*nVertex);
hist(tmpStats, -80:10:80)

xlabel('Optimization Best Error - Optimization Error');
ylabel('Frequency')
