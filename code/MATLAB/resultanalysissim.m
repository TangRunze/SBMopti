close all
clear all

%% --- Simulation ---

r = -1;

% pathExtra = '/';
pathExtra = '/latest/';

% Parameter Setting
nVertex = 150;
gStart = 1;
gEnd = 200;

rho = [1/3, 1/3, 1/3];


B = [0.4, 0.2, 0.2; 0.2, 0.4, 0.2; 0.2, 0.2, 0.4];

lambdaVec = [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999];
% lambdaVec = [0.5];
nLambda = length(lambdaVec);

nBlock = 3;
dimLatentPosition = 3;
isGMM = 0;



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
errorOptiOracleVec = zeros(1, maxIter);
errorASGEGMMVec = zeros(1, maxIter);
loglikASGEVec = zeros(1, maxIter);
loglikOptiVec = zeros(1, maxIter);
loglikOptiOracleVec = zeros(1, maxIter);


tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    100, 'MaxFunEvals', 5000);



for iInd = 1:maxIter
    iGraph = ind(iInd);
    errorRateASGETmp = [];
    errorRateOptiTmp = [];
    loglikASGETmp = [];
    loglikOptiTmp = [];
    for iLambda = 1:nLambda
        lambda = lambdaVec(iLambda);
        if exist(['./results' pathExtra 'results-SBMopti-Block-sim-dir-n' ...
                num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
                num2str(B(1, 2)) '-GMM' num2str(isGMM) '-r' num2str(r) ...
                '-lambda' num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
            load(['./results' pathExtra 'results-SBMopti-Block-sim-dir-n' ...
                num2str(nVertex) '-diag' num2str(B(1, 1)) '-offdiag' ...
                num2str(B(1, 2)) '-GMM' num2str(isGMM) '-r' num2str(r) ...
                '-lambda' num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
            errorRateASGETmp = [errorRateASGETmp errorRateASGE];
            errorRateOptiTmp = [errorRateOptiTmp errorRateOpti];
            loglikASGETmp = [loglikASGETmp loglikASGE];
            loglikOptiTmp = [loglikOptiTmp loglikOpti];
        end
    end
    [maxLoglikOpti, indLoglikOpti] = max(loglikOptiTmp);
    [minErrorRateOpti, indErrorRateOpti] = min(errorRateOptiTmp);
    errorASGEVec(iInd) = mean(errorRateASGE);
    errorOptiVec(iInd) = errorRateOptiTmp(indLoglikOpti);
    errorOptiOracleVec(iInd) = minErrorRateOpti;
    loglikASGEVec(iInd) = mean(loglikASGE);
    loglikOptiVec(iInd) = maxLoglikOpti;
    loglikOptiOracleVec(iInd) = loglikOptiTmp(indErrorRateOpti);
    
    
    [adjMatrix, adjMatrixDA] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, ...
        rho, tauStar, r, iGraph, projectoptions);
    
    xHat0 = asge(adjMatrixDA, dimLatentPosition);
    %     xHat0 = asge(adjMatrix, dimLatentPosition);
    
    [tauHat0, nuHat0] = clusterX(xHat0, nBlock, 1);
    errorASGEGMMVec(iInd) = errorratecalculator(tauStar, tauHat0, nVertex, nBlock);
    
end


errorMeanASGEGMM = mean(errorASGEGMMVec);
errorCIASGEGMM = [errorMeanASGEGMM - 1.96*std(errorASGEGMMVec)/sqrt(maxIter), ...
    errorMeanASGEGMM + 1.96*std(errorASGEGMMVec)/sqrt(maxIter)];
errorMedianASGEGMM = median(errorASGEGMMVec);


errorMeanASGE = mean(errorASGEVec);
errorCIASGE = [errorMeanASGE - 1.96*std(errorASGEVec)/sqrt(maxIter), ...
    errorMeanASGE + 1.96*std(errorASGEVec)/sqrt(maxIter)];
errorMedianASGE = median(errorASGEVec);

errorMeanOpti = mean(errorOptiVec);
errorCIOpti = [errorMeanOpti - 1.96*std(errorOptiVec)/sqrt(maxIter), ...
    errorMeanOpti + 1.96*std(errorOptiVec)/sqrt(maxIter)];
errorMedianOpti = median(errorOptiVec);

errorMeanOptiOracle = mean(errorOptiOracleVec);
errorCIOptiOracle = [errorMeanOptiOracle - 1.96*std(errorOptiOracleVec)/sqrt(maxIter), ...
    errorMeanOptiOracle + 1.96*std(errorOptiOracleVec)/sqrt(maxIter)];
errorMedianOptiOracle = median(errorOptiOracleVec);



[errorMeanASGEGMM, errorMeanASGE, errorMeanOpti, errorMeanOptiOracle]

[errorMedianASGEGMM, errorMedianASGE, errorMedianOpti, errorMedianOptiOracle]

% [loglikMeanASGE, loglikMeanOpti]

%% --- Paired Test ---

% Alternative hypothesis: Opti < ASGE.

errorPair = zeros(2, maxIter);
errorPair(1, :) = errorOptiVec;
errorPair(2, :) = errorASGEVec;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) < errorPair(2, :));
pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) == errorPair(2, :)), 0.5)


% Alternative hypothesis: OptiOracle < ASGE.

errorPair = zeros(2, maxIter);
errorPair(1, :) = errorOptiOracleVec;
errorPair(2, :) = errorASGEVec;

% 1-sided sign-test
tmpStats = sum(errorPair(1, :) < errorPair(2, :));
pValue = 1 - binocdf(tmpStats - 1, maxIter - ...
    sum(errorPair(1, :) == errorPair(2, :)), 0.5)



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





%% --- Plot ---

ind = [1:8];
Perr =  [errorRateASGEVecMean(ind)';
    errorRateOptiVecMean(ind)';
    errorRateOptiOptiVecMean(ind)'];

%flat
lower_opti = errorRateOptiVecMean(ind)' - errorRateOptiVecVar(ind)';
upper_opti = errorRateOptiVecMean(ind)' + errorRateOptiVecVar(ind)';
CI_opti = [lower_opti;upper_opti];
%asge
lower_asge = errorRateASGEVecMean(ind)' - errorRateASGEVecVar(ind)';
upper_asge = errorRateASGEVecMean(ind)' + errorRateASGEVecVar(ind)';
CI_asge = [lower_asge;upper_asge];
%Oracle Bayes
lower_oracle_opti = errorRateOptiOptiVecMean(ind)' - errorRateOptiVecVar(ind)';
upper_oracle_opti = errorRateOptiOptiVecMean(ind)' + errorRateOptiOptiVecVar(ind)';
CI_oracle_opti = [lower_oracle_opti;upper_oracle_opti];

% normal scale

figure
set(gcf,'Color',[1,1,1])
plot(Perr(1,:),'b-','LineWidth',2);
axis([0 8 0 0.7]);
hold on

plot(Perr(2,:),'r-','LineWidth',2);
plot(Perr(3,:),'g-','LineWidth',2);
hold all
leg = legend('ASGE','Opti','Oracle Opti');
legend('boxoff')

set(leg,'FontSize',15)
set(gca,'ygrid','on')
set(gca,'Xtick',1:98,'XTickLabel',{'0', '1', '5', '10', '20', '50', '100', 'Infinity'})

plotshaded(1:size(Perr, 2), CI_asge,'b')
plotshaded(1:size(Perr, 2), CI_opti,'r')
plotshaded(1:size(Perr, 2), CI_oracle_opti,'g')
hold off
xlabel('r - Paramter of Dirichlet','FontSize',12), ylabel('Classification Error','FontSize',12)

% print(gcf,'-dpng','-r600','fig1.png')





