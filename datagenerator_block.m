function [pMatrix, muHat, sigmaHat, tauHat, pTauHat] = ...
    datagenerator_block(nVertex, nBlock, dimLatentPosition, B, rho, ...
    tauStar, epsilonInB, r, iProbMatrix, projectoptions)
% Generate data if there does not exist one, otherwise read the
% existing data.

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
        xDirichlet = [xDirichlet; drchrnd(r*nuStar(iBlock, :), nVec(iBlock))];
    end
    
    % Get the probability matrix with Dirchlet samples.
    pMatrix = xDirichlet*xDirichlet';
    
    % Obtain estimates from ASGE o GMM.
    xHat = asge(pMatrix, dimLatentPosition);
    
    gm = fitgmdist(xHat, nBlock, 'Replicates', 10, 'Regularize', 1e-12);
    
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
    save(['data/sim-dir-n' num2str(nVertex) '-eps' num2str(epsilonInB) ...
        '-r' num2str(r) '-pmatrix' int2str(iProbMatrix) '.mat'],...
        'pMatrix', 'tauHat', 'pTauHat', 'muHat', 'sigmaHat');
else
    % Read the existing data
    data = load(['data/sim-dir-n' num2str(nVertex) '-eps' ...
        num2str(epsilonInB) '-r' num2str(r) '-pmatrix' ...
        int2str(iProbMatrix) '.mat']);
    pMatrix = data.pMatrix;
    muHat = data.muHat;
    sigmaHat = data.sigmaHat;
    tauHat = data.tauHat;
    pTauHat = data.pTauHat;
end
