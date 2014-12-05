close all
clear all

%% --- Parameters Setting ---

nVertex = 150;
nBlock = 3;
dimLatentPosition = nBlock;
diag = 0.5;
offdiag = 0.1;
r = 5;
gStart = 1;
gEnd = 1;

lambdaVec = (1:999)/1000;
% lambdaVec = [0.1];
nLambda = length(lambdaVec);

%% --- Result of Error Rates ---

nuVec = zeros(nLambda, nBlock, dimLatentPosition);
deltaNu = zeros(1, nLambda - 1);
resultVec = [];

% Go over different lambda.
for iLambda = 1:nLambda
    lambda = lambdaVec(iLambda);
    iProbMatrix = gStart;
    if exist(['./results/results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
            num2str(offdiag) '-r' num2str(r) '-lambda' ...
            num2str(lambda) '-pmatrix' num2str(iProbMatrix) '.mat'])
        load(['./results/results-SBMopti-Block-sim-dir-n' ...
            num2str(nVertex) '-diag' num2str(diag) '-offdiag' ...
            num2str(offdiag) '-r' num2str(r) '-lambda' ...
            num2str(lambda) '-pmatrix' num2str(iProbMatrix) '.mat']);
        nuVec(iLambda, :, :) = nuBest;
        resultVec = [resultVec iLambda];
    end
    if (iLambda > 1)
        deltaNu(iLambda - 1) = norm(squeeze(nuVec(iLambda, :, :)) - ...
            squeeze(nuVec(iLambda - 1, :, :)));
    end
end


nVec = 1:100;

% plot(deltaNu);
plot3(squeeze(nuVec(:, 1, 1)), squeeze(nuVec(:, 1, 2)), squeeze(nuVec(:, 1, 3)), '.b');
hold on;
plot3(squeeze(nuVec(:, 2, 1)), squeeze(nuVec(:, 2, 2)), squeeze(nuVec(:, 2, 3)), '.b');
hold on;
plot3(squeeze(nuVec(:, 3, 1)), squeeze(nuVec(:, 3, 2)), squeeze(nuVec(:, 3, 3)), '.b');
hold on;
plot3(squeeze(nuVec(nVec, 1, 1)), squeeze(nuVec(nVec, 1, 2)), squeeze(nuVec(nVec, 1, 3)), '.r');
hold on;
plot3(squeeze(nuVec(nVec, 2, 1)), squeeze(nuVec(nVec, 2, 2)), squeeze(nuVec(nVec, 2, 3)), '.r');
hold on;
plot3(squeeze(nuVec(nVec, 3, 1)), squeeze(nuVec(nVec, 3, 2)), squeeze(nuVec(nVec, 3, 3)), '.r');
hold on;


% nuStar

nBlock = 3;
dimLatentPosition = nBlock;
epsilonInB = 0.2;
rho = repmat(1/nBlock, 1, nBlock);
tol = 1e-4;
maxIter = 100;
muB = 0.3;
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);


nuStar = asge(B, dimLatentPosition);

% Pre-projection
rotationMatrix = fmincon(@(x) projectobjectivefun(x, ...
    dimLatentPosition, nuStar), ...
    reshape(eye(dimLatentPosition), dimLatentPosition^2, 1), ...
    [], [], [], [], - ones(dimLatentPosition^2, 1), ...
    ones(dimLatentPosition^2, 1), @(x) ...
    projectconditionfun(x, dimLatentPosition, dimLatentPosition, ...
    nuStar), projectoptions);
rotationMatrix = reshape(rotationMatrix, dimLatentPosition, ...
    dimLatentPosition);

% Rotate the latent positions
nuStar = nuStar*rotationMatrix;


hold on;
plot3(nuStar(1, 1), nuStar(1, 2), nuStar(1, 3), 'go')
plot3(nuStar(2, 1), nuStar(2, 2), nuStar(2, 3), 'go')
plot3(nuStar(3, 1), nuStar(3, 2), nuStar(3, 3), 'go')




% squeeze(nuVec(resultVec(1), 1, :))'
% nuStar
% (r*nuStar + 1)/(r+4)