close all
clear all

nVertex = 150;
nBlock = 3;
dimLatentPosition = nBlock;
muB = 0.3;
epsilonInB = 0.2;
diag = muB + epsilonInB;
offdiag = muB - epsilonInB;
gStart = 1;
gEnd = 500;

rho = repmat(1/nBlock, 1, nBlock);
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

rVec = [1 3 5 6 7 8 9 10 15 20];
nR = length(rVec);

lambdaVec = [0.001 0.01 0.1 0.9 0.99 0.999];
nLambda = length(lambdaVec);



% Go over different r.
for iR = 1:nR
    r = rVec(iR);
    % Go over all Monte Carlo replicates.
    for iGraph = gStart:gEnd
        % Go over different lambda.
        for iLambda = 1:nLambda
            lambda = lambdaVec(iLambda);
            if exist(['./results/results-SBMopti-Block-sim-dir-n' ...
                    num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                    num2str(offdiag) '-r' num2str(r) '-lambda' ...
                    num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat'])
                load(['./results/results-SBMopti-Block-sim-dir-n' ...
                    num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
                    num2str(offdiag) '-r' num2str(r) '-lambda' ...
                    num2str(lambda) '-adjMatrix' num2str(iGraph) '.mat']);
                csvwrite(['results/KLdata/xBest-n' num2str(nVertex) ...
                    '-diag' num2str(diag) '-offdiag' num2str(offdiag) ...
                    '-r' num2str(r) '-lambda' num2str(lambda) ...
                    '-adjMatrix' num2str(iGraph)], xBest);
            end
        end
    end
end





