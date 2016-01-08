function plot2Dxhat(xHat, muHat, tauStar, nuStar, muHat2)

nv1 = (tauStar == 1);
nv2 = (tauStar == 2);

% axis([-1 1 -1 1 -1 1])

plot(xHat(nv1, 1), xHat(nv1, 2), 'g.','markersize',12);
hold on;
plot(xHat(nv2, 1), xHat(nv2, 2), 'b.','markersize',12);

plot(muHat(:, 1), muHat(:, 2), 'ro','markersize',12,'linewidth',8);

if nargin>3
    plot(nuStar(:, 1), nuStar(:, 2), 'k+','markersize',18,'linewidth',8);
end

if nargin>4
    plot(muHat2(:, 1), muHat2(:, 2), 'yx','markersize',18,'linewidth',8);
end

legend('1','2','k','truth','g')