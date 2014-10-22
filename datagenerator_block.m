function [pMatrix, muHat, sigmaHat, tauHat, pTauHat] = ...
    datagenerator_block(nVertex, nBlock, dimLatentPosition, B, rho, ...
    tauStar, epsilonInB, delta, iProbMatrix)
% Generate data if there does not exist one, otherwise read the
% existing data.

if exist(['data/sim-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
        '-delta' num2str(delta) '-pmatrix' int2str(iProbMatrix) '.mat'], ...
        'file') == 0
    
    disp(['Generating Probability Matrix ' int2str(iProbMatrix) '...'])
    
    % Get the point-mass latent positions of a pure block model.
    nuStar = asge(B, dimLatentPosition);
    
    % Generate the noise corresponding to each vertex.
    noiseTmp = unifrnd(-delta, delta, 2*nVertex, dimLatentPosition);
    nv = (sum(noiseTmp.^2, 2) <= delta^2);
    noiseVector = noiseTmp(nv, :);
    while (size(noiseVector, 1) < nVertex)
        noiseTmp = unifrnd(0, delta, 2*nVertex, dimLatentPosition);
        nv = (sum(noiseTmp.^2, 2) <= delta^2);
        noiseVector = [noiseVector; noiseTmp(nv, :)];
    end
    noiseVector = noiseVector(1:nVertex, :);
    
    % Get the probability matrix with noise delta.
    xNoise = nuStar(tauStar, :) + noiseVector;
    pMatrix = xNoise*xNoise';
    
    % Obtain estimates from ASGE o GMM.
    xHat = asge(pMatrix, dimLatentPosition);
    
    gm = fitgmdist(xHat, nBlock, 'Replicates', 10);
    
    tauHat = cluster(gm, xHat)';
    % pihat = gm.PComponents;
    pTauHat = posterior(gm, xHat)';
    muHat = gm.mu;
    sigmaHat = gm.Sigma;
    
    % Plot
    % cl_nv = false(K,n);
    % for i = 1:K
    %     cl_nv(i,:) = (idx == i);
    % end
    % scatter(Xhat(cl_nv(1,:),1),Xhat(cl_nv(1,:),2),10,'r+');
    % hold on
    % scatter(Xhat(cl_nv(2,:),1),Xhat(cl_nv(2,:),2),10,'bo');
    % hold on
    % scatter(Xhat(cl_nv(3,:),1),Xhat(cl_nv(3,:),2),10,'g.');
    % hold off
    % legend('Cluster 1','Cluster 2','Cluster 3','Location','NW')
    
    % Save the data
    save(['data/sim-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
        '-delta' num2str(delta) '-pmatrix' int2str(iProbMatrix) '.mat'],...
        'pMatrix', 'tauHat', 'pTauHat', 'muHat', 'sigmaHat');
else
    % Read the existing data
    data = load(['data/sim-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
        '-delta' num2str(delta) '-pmatrix' int2str(iProbMatrix) '.mat']);
    pMatrix = data.pMatrix;
    muHat = data.muHat;
    sigmaHat = data.sigmaHat;
    tauHat = data.tauHat;
    pTauHat = data.pTauHat;
end
