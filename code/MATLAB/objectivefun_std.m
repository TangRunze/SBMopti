function val = objectivefun_std(x, adjMatrix, nuHat, lambda, ...
    nVertex, dimLatentPosition, tauHat)
% Objective functions in optimization problem

X = reshape(x, nVertex, dimLatentPosition);

% val = (1 - lambda) * norm(adjMatrix - X*X', 'fro')^2 + ...
%     lambda*norm(X - nuHat(tauHat, :), 'fro')^2 / dimLatentPosition;

adjMatrix(1:(nVertex+1):end) = diag(X*X');

val = (1 - lambda) * norm(adjMatrix - X*X', 'fro')^2/(nVertex*(nVertex-1)) + ...
    lambda*norm(X - nuHat(tauHat, :), 'fro')^2/(nVertex*dimLatentPosition);

% val = (1 - lambda) * norm(adjMatrix - X*X', 'fro')^2 + ...
%     lambda*norm(X - nuHat(tauHat, :), 'fro')^2;

% grad = - 4*(adjMatrix + X*X')*X + 2*lambda*(X + nuHat(tauHat, :));

% val = - 2*trace(adjMatrix*(X*X')) + trace((X'*X)*(X'*X)) + lambda*trace((X - nuHat(tauHat, :))'*(X - nuHat(tauHat, :)));

end