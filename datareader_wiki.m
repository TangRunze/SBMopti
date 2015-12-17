function [nVertex, adjMatrix, adjMatrixDA, ind, tauStar] = ...
    datareader_wiki(iGraph)
% Read the existing data.

rng shuffle;

if exist(['data/real-graph' int2str(iGraph) '.mat'], 'file') == 0
    
    disp(['Generating subgraph ' int2str(iGraph) '...'])
    
    % Read the entire graph
    edgelist = textread(['data/edgelist']);
    nVertex0 = max(max(edgelist));
    adjMatrix0 = zeros(nVertex0);
    for iEdgelist = 1:size(edgelist, 1)
        adjMatrix0(edgelist(iEdgelist, 1), edgelist(iEdgelist, 2)) = 1;
        adjMatrix0(edgelist(iEdgelist, 2), edgelist(iEdgelist, 1)) = 1;
    end
    
    % Read the labels
    tauStar0 = textread(['data/classes'])';
    
    % Make it symmetric
    adjMatrix0 = adjMatrix0 + adjMatrix0';
    adjMatrix0 = 1*(adjMatrix0 > 0);
    
    t1 = 1;
    t2 = 2;
    t3 = 3;
    
    % Take the first 3 classes
    nv123 = (tauStar0 == t1) | (tauStar0 == t2) | (tauStar0 == t3);
    adjMatrix123 = adjMatrix0(nv123, nv123);
    tauStar123 = tauStar0(nv123);
    
    tauTmp = tauStar123;
    tauStar123(tauTmp == t1) = 1;
    tauStar123(tauTmp == t2) = 2;
    tauStar123(tauTmp == t3) = 3;
    
    % Excluding the singletons
    [~, C] = graphconncomp(sparse(adjMatrix123));
    tmp = tabulate(C);
    [~, indexLCC] = max(tmp(:, 3));
    nvLCC = (C == indexLCC);
    
    adjMatrix123 = adjMatrix123(nvLCC, nvLCC);
    tauStar123 = tauStar123(nvLCC);
    nVertex123 = length(tauStar123);
    
    % Sort by classes
    nv1 = (tauStar123 == 1);
    nv2 = (tauStar123 == 2);
    nv3 = (tauStar123 == 3);
    n1 = sum(nv1);
    n2 = sum(nv2);
    n3 = sum(nv3);
    index123 = 1:nVertex123;
    adjMatrix123 = adjMatrix123([index123(nv1), index123(nv2), ...
        index123(nv3)], [index123(nv1), index123(nv2), index123(nv3)]);
    tauStar123 = tauStar123([index123(nv1), index123(nv2), index123(nv3)]);
    
    % Diagonal Augmentation.
    adjMatrixDA = adjMatrix123 + diag(sum(adjMatrix123))/...
        (size(adjMatrix123, 1) - 1);
    
    % Take the subgraph randomly
    ind1 = randsample(n1, 100, true);
    ind2 = randsample(n2, 100, true) + n1;
    ind3 = randsample(n3, 100, true) + n1 + n2;
    ind = [ind1; ind2; ind3];
    
    adjMatrixDA = adjMatrixDA(ind, ind);
    
    adjMatrix = adjMatrix123(ind, ind);
    tauStar = tauStar123(ind);
    nVertex = length(ind);
    
    % Save the data
    save(['data/real-graph' int2str(iGraph) '.mat'], 'adjMatrixDA', ...
        'ind', 'adjMatrix', 'nVertex', 'tauStar');
    
else
    % Read the existing data
    data = load(['data/real-graph' int2str(iGraph) '.mat']);
    adjMatrixDA = data.adjMatrixDA;
    ind = data.ind;
    adjMatrix = data.adjMatrix;
    nVertex = data.nVertex;
    tauStar = data.tauStar;
end
