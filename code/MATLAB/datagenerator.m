function [adjMatrix, adjMatrixDA,nuStar] = ...
    datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, ...
    tauStar, r, iGraph, projectoptions)
% Generate data if there does not exist one, otherwise read the
% existing data.

if exist(['data/sim-dir-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-r' num2str(r) '-graph' ...
        int2str(iGraph) '.mat'], 'file') == 0
    
    disp(['Generating Adjacency Matrix ' int2str(iGraph) '...'])
    
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
    % latent positions (actually with bias).
    if (r ~= -1)
        nVec = nVertex * rho;
        xDirichlet = [];
        for iBlock = 1:nBlock
            xTmp = drchrnd(r*[nuStar(iBlock, :), ...
                1 - sum(nuStar(iBlock, :))] + ...
                ones(1, dimLatentPosition + 1), nVec(iBlock));
            xDirichlet = [xDirichlet; xTmp(:, 1:dimLatentPosition)];
        end
    else
        xDirichlet = nuStar(tauStar, :);
    end
    
    % Get the probability matrix with Dirchlet samples.
    pMatrix = xDirichlet*xDirichlet';
    
    adjMatrix = reshape(binornd(ones(1, nVertex*nVertex), ...
        reshape(pMatrix, 1, nVertex*nVertex)), nVertex, nVertex);
    adjMatrix = triu(adjMatrix, 1);
    adjMatrix = adjMatrix + adjMatrix';
    
    % Diagonal Augmentation.
    adjMatrixDA = adjMatrix + diag(sum(adjMatrix))/...
        (size(adjMatrix, 1) - 1);
    
    % Save the data
    save(['../../data/sim-dir-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-r' num2str(r) '-graph' ...
        int2str(iGraph) '.mat'], 'adjMatrixDA', 'adjMatrix', 'nuStar');
else
    % Read the existing data
    data = load(['../../data/sim-dir-n' num2str(nVertex) '-diag' num2str(B(1, 1)) ...
        '-offdiag' num2str(B(1, 2)) '-r' num2str(r) '-graph' ...
        int2str(iGraph) '.mat']);
    adjMatrixDA = data.adjMatrixDA;
    adjMatrix = data.adjMatrix;
    nuStar = data.nuStar;
end
