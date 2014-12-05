
%% --- Quick Setting ---
nVertex = 150;
nBlock = 3;
dimLatentPosition = nBlock;
epsilonInB = 0.2;
rho = repmat(1/nBlock, 1, nBlock);
tol = 1e-4;
maxIter = 100;
muB = 0.5;
r = 100;
B = (muB - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end
options = optimoptions('fmincon', 'TolX', 1e-6, ...
    'MaxIter', 10000, 'MaxFunEvals', 10000, ...
    'Algorithm', 'interior-point');
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);

%% --- Parameter Setting ---
% nVertex selects the number of vertices in the graph.
% nVertex = 150;

% nBlock selects the number of blocks in the stochastic blockmodel.
% nBlock = 3;

% dimLatentPosition selects the dimension of latent positions.
dimLatentPosition = nBlock;

% true block proportion
rho = repmat(1/nBlock, 1, nBlock);

% epsilonInB controls the true model. The probability matrix
%       B = (0.5 - epsilonInB)*J + 2*epsilonInB*I
% epsilonInB should be inside [0, 0.5].
% epsilonInB = 0.1;

% block probability matrix
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);

% true tau_star (1-by-nVertex)
tauStar = [];
nVectorStar = nVertex*rho;
for i = 1:nBlock
    tauStar = [tauStar, i*ones(1, nVectorStar(i))];
end

%% --- Optimization Setting ---
options = optimoptions('fmincon', 'TolX', 1e-6, ...
    'MaxIter', 10000, 'MaxFunEvals', 10000, ...
    'Algorithm', 'interior-point'); %, 'GradObj', 'on');
% Algorithm: 'interior-point', 'trust-region-reflective', 'sqp', 'active-set'
projectoptions = optimoptions('fmincon', 'TolX', 1e-6, 'MaxIter', ...
    10000, 'MaxFunEvals', 10000);

lambdaVec = [0.001 0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.999];

for iProbMatrix = gStart:gEnd
    
    if exist(['data/sim-dir-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
            '-r' num2str(r) '-pmatrix' int2str(iProbMatrix) '.mat'], ...
            'file') == 0
        
        disp(['Generating Probability Matrix ' int2str(iProbMatrix) '...'])
        
        % Get the point-mass latent positions of a pure block model.
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
        
        % Generate samples from a Dirichlet distribution centered at true
        % latent positions.
        nVec = nVertex * rho;
        xDirichlet = [];
        for iBlock = 1:nBlock
            xDirichlet = [xDirichlet; ...
                drchrnd(r*nuStar(iBlock, :) + ones(1, dimLatentPosition), nVec(iBlock))];
        end
        
        % Get the probability matrix with Dirchlet samples.
        pMatrix = xDirichlet*xDirichlet';
        
        % Obtain estimates from ASGE o GMM.
        [V, D] = eigs(pMatrix, dimLatentPosition, 'LA');
        xHat = V*sqrt(D);
        
        
        gm = fitgmdist(xHat, nBlock, 'Replicates', 10, 'Regularize', 1e-12);
        
        tauHat = cluster(gm, xHat)';
        % pihat = gm.PComponents;
        pTauHat = posterior(gm, xHat)';
        muHat = gm.mu;
        sigmaHat = gm.Sigma;
    end
    
    
end

%% --- Close Parallel Computing ---
delete(gcp('nocreate'))

