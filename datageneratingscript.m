
gStart = 1;
gEnd = 1000;
nVertex = 150;
nBlock = 3;
dimLatentPosition = 3;
epsilonInB = 0.1;
B = (0.5 - epsilonInB)*ones(nBlock, nBlock) + 2*epsilonInB*eye(nBlock);
rho = repmat(1/nBlock, 1, nBlock);


% delete(gcp('nocreate'))
% parpool(12);

for iGraph = gStart:gEnd
    
    %% --- Generate/Read Data ---
    % Generate data if there does not exist one, otherwise read the
    % existing data.
    [adjMatrix, nuHat, ~, tauHat, ~] = ...
        datagenerator(nVertex, nBlock, dimLatentPosition, B, rho, ...
        epsilonInB, iGraph);
end

% delete(gcp('nocreate'))
