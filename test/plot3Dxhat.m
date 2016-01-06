function plot3Dxhat(xHat, muHat, tauStar)

nv1 = (tauStar == 1);
nv2 = (tauStar == 2);
nv3 = (tauStar == 3);

figure;
axis([-1 1 -1 1 -1 1])

plot3(xHat(nv1, 1), xHat(nv1, 2), xHat(nv1, 3), 'g.');
hold on;
plot3(xHat(nv2, 1), xHat(nv2, 2), xHat(nv2, 3), 'b.');
plot3(xHat(nv3, 1), xHat(nv3, 2), xHat(nv3, 3), 'm.');

plot3(muHat(:, 1), muHat(:, 2), muHat(:, 3), 'ro');