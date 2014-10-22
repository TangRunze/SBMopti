function val = objectivefun_std_lambda1(adjMatrix, nuHat, ...
    dimLatentPosition, nBlock, tauHat)
% Objective functions in optimization problem

nu = reshape(nuHat, nBlock, dimLatentPosition);
X = nu(tauHat, :);

val = norm(adjMatrix - X*X', 'fro')^2;

% grad = - 4*(adjMatrix + X*X')*X + 2*lambda*(X + nuHat(tauHat, :));

% val = - 2*trace(adjMatrix*(X*X')) + trace((X'*X)*(X'*X)) + lambda*trace((X - nuHat(tauHat, :))'*(X - nuHat(tauHat, :)));

end