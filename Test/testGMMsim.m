clear all
close all

nGraph = 500;
nVertex = 150;
eps = 0.1;
nBlock = 3;
dimLatentPosition = nBlock;

rho = repmat(1/nBlock, 1, nBlock);
B = (0.5 - eps)*ones(nBlock) + 2*eps*eye(nBlock);
tauStar = [ones(1,nVertex/nBlock), 2*ones(1,nVertex/nBlock), 3*ones(1,nVertex/nBlock)];

nVectorStar = nVertex*rho;
nVectorStarStart = cumsum(nVectorStar);
nVectorStarStart = [1, nVectorStarStart(1:(end-1)) + 1];
nVectorStarEnd = cumsum(nVectorStar);

for iGraph = 1:nGraph
    
    disp(['Generating graph ' int2str(iGraph) '...'])
    
    % Part 1: Generate graph adjacency matrix.
    adjMatrix = zeros(nVertex);
    for iBlock = 1:nBlock
        for jBlock = iBlock:nBlock
            adjMatrix(nVectorStarStart(iBlock):nVectorStarEnd(iBlock), ...
                nVectorStarStart(jBlock):nVectorStarEnd(jBlock)) = ...
                binornd(1, B(iBlock,jBlock), nVectorStar(iBlock), ...
                nVectorStar(jBlock));
        end
    end
    adjMatrix = triu(adjMatrix, 1);
    adjMatrix = adjMatrix + adjMatrix';
    
    fileID = fopen(['./testdata/sim-graph', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',adjMatrix(:));
    fclose(fileID);
    
    xHat = asge(adjMatrix, dimLatentPosition);
    
    fileID = fopen(['./testdata/sim-xhat', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',xHat(:));
    fclose(fileID);
end


errorGMM = zeros(1, nGraph);

for iGraph = 1:nGraph
    % Part 2: Obtain estimates from ASGE o GMM.
    
    disp(['Generating graph ' int2str(iGraph) '...'])
    
    tauInit = textread(['./testdata/sim-tauinit', num2str(iGraph), '-n',...
        num2str(nVertex), '-eps', num2str(eps)]);
    
    xHat = textread(['./testdata/sim-xhat', num2str(iGraph), '-n',...
        num2str(nVertex), '-eps', num2str(eps)]);
    xHat = reshape(xHat, nVertex, dimLatentPosition);
    
    gm = fitgmdist(xHat, nBlock, 'Start', tauInit);
    
    % gm = fitgmdist(xHat, nBlock, 'Start', tauInit, 'Options', ...
    %     statset('MaxIter', 350, 'TolFun', sqrt(2.22e-16)*100));
    % gm = fitgmdist(xHat, nBlock, 'Start', tauInit, 'Options', ...
    %     statset('MaxIter', 350));
    % gm = fitgmdist(xHat, nBlock, 'Replicates', 10);
    
    tauHat = cluster(gm, xHat)';
    % pihat = gm.PComponents;
    pTauHat = posterior(gm, xHat)';
    muHat = gm.mu;
    sigmaHat = gm.Sigma;
    proHat = gm.PComponents;
    
    fileID = fopen(['./testdata/sim-tauhat', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',tauHat(:));
    fclose(fileID);
    
    fileID = fopen(['./testdata/sim-ptauhat', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',pTauHat(:));
    fclose(fileID);
    
    fileID = fopen(['./testdata/sim-muhat', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',muHat(:));
    fclose(fileID);
    
    fileID = fopen(['./testdata/sim-sigmahat', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',sigmaHat(:));
    fclose(fileID);
    
    fileID = fopen(['./testdata/sim-prohat', num2str(iGraph), '-n', ...
        num2str(nVertex), '-eps', num2str(eps)],'w');
    fprintf(fileID,'%d\n',proHat(:));
    fclose(fileID);
    
    errorGMM(iGraph) = errorratecalculator(tauStar, tauHat, nVertex, nBlock);

end