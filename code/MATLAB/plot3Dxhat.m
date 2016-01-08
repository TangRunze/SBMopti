function plot3Dxhat(xHat, muHat, tauStar, muHat2)

nv1 = (tauStar == 1);
nv2 = (tauStar == 2);
nv3 = (tauStar == 3);

axis([-1 1 -1 1 -1 1])

plot3(xHat(nv1, 1), xHat(nv1, 2), xHat(nv1, 3), 'g.','markersize',12);
hold on;
plot3(xHat(nv2, 1), xHat(nv2, 2), xHat(nv2, 3), 'b.','markersize',12);
plot3(xHat(nv3, 1), xHat(nv3, 2), xHat(nv3, 3), 'm.','markersize',12);

plot3(muHat(:, 1), muHat(:, 2), muHat(:, 3), 'ro','markersize',12,'linewidth',8);

if nargin==4
plot3(muHat2(:, 1), muHat2(:, 2), muHat2(:, 3), 'k+','markersize',18,'linewidth',8);
end