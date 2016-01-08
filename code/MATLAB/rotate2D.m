function [Z, R] = rotate2D(X, Y)

[U, ~, V] = svd(X'*Y);
R = U*V';
Z = X*R;
