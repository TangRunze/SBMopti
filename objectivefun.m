function [val, grad] = objectivefun(x, adjMatrix, nuHat, lambda, ...
    nVertex, dimLatentPosition, tauHat)
% Objective functions in optimization problem

X = reshape(x, nVertex, dimLatentPosition);

val = norm(adjMatrix - X*X', 'fro')^2 + ...
    lambda*norm(X - nuHat(tauHat, :), 'fro')^2;

grad = - 4*(adjMatrix + X*X')*X + 2*lambda*(X + nuHat(tauHat, :));

% val = - 2*trace(adjMatrix*(X*X')) + trace((X'*X)*(X'*X)) + lambda*trace((X - nuHat(tauHat, :))'*(X - nuHat(tauHat, :)));

end