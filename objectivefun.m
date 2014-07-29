function val = objectivefun(x, adjMatrix, nuHat, lambda, nVertex, ...
    dimLatentPosition, tauHat)
% Objective functions in optimization problem

X = reshape(x, nVertex, dimLatentPosition);

% penalty = 0;
% for iVertex = 1:nVertex
%     penalty = penalty + norm(X(iVertex, :) - nuHat(tauHat(iVertex), :), 2);
% end
% penalty = penalty*lambda;
% 
% val = norm(adjMatrix - X*X', 'fro') + penalty;

val = norm(adjMatrix - X*X', 'fro')^2 + ...
    lambda*norm(X - nuHat(tauHat, :), 'fro')^2;

% val = norm(adjMatrix - X*X', 'fro');

end