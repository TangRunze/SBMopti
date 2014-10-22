rng('default')
n = 1000;
X = mvnrnd([0.4, 0.4], 0.01*eye(2), n);

theta = pi/2;
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
Y = X*R;

X0 = 0.4*ones(n, 2);

[d,Z,tr] = procrustes(X0, Y, 'scaling', false, 'reflection', false);

norm(Z - (tr.b*Y*tr.T + tr.c))

plot(X(:,1),X(:,2),'rx',Y(:,1),Y(:,2),'b.',Z(:,1),Z(:,2),'bx');
% plot(X(:, 1), X(:, 2), 'rx', Y(:,1), Y(:,2), 'b.');



% plot(X(:,1),X(:,2),'rx')