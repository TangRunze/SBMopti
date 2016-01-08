function plot2Dxhatheatmap(xHat, tauStar, nuStar, nuHat_k, nuHat_g, ...
    sigmaHat_g, proportionHat_g)

%% Rotate nuStar
nuStar = rotate2D(nuStar, nuHat_g);


%% Heat Map
gridSize = 0.01;
[x, y] = meshgrid(-1:gridSize:1, -1:gridSize:1);

l = size(x, 1);
pointTmp = [reshape(x, l^2, 1), reshape(y, l^2, 1)];

pdfGaussian = gmmpdfcalculator(pointTmp, nuHat_g, sigmaHat_g, ...
    proportionHat_g);

z = reshape(pdfGaussian, l, l);

imagesc([-1, 1], [-1, 1], z, 'CDataMapping', 'scaled');
hold on;
colorbar;


%% Plot other Features
nv1 = (tauStar == 1);
nv2 = (tauStar == 2);

markerSize = 5;

plot(xHat(nv1, 1), xHat(nv1, 2), 'g.', 'markersize', markerSize);
plot(xHat(nv2, 1), xHat(nv2, 2), 'w.', 'markersize', markerSize);

plot(nuHat_k(:, 1), nuHat_k(:, 2), 'ro', 'markersize', markerSize, ...
    'linewidth', markerSize);

plot(nuHat_g(:, 1), nuHat_g(:, 2), 'mx', 'markersize', markerSize, ...
    'linewidth', markerSize);

plot(nuStar(:, 1), nuStar(:, 2), 'k+', 'markersize', markerSize, ...
    'linewidth', markerSize);

legend('1', '2', 'k', 'g', 'truth')


%% Save the Figure
saveas(gca, '../../results/2DHeatMap.eps', 'epsc');
