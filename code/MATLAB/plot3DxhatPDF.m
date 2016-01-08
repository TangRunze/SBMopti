function plot3DxhatPDF(xHat, tauStar, nuStar, nuHat_k, nuHat_g, ...
    sigmaHat_g, proportionHat_g)


nuStar = rotate2D(nuStar, nuHat_g);

%% Plot Components

nv1 = (tauStar == 1);
nv2 = (tauStar == 2);


hold on;
plot3(xHat(nv1, 1), xHat(nv1, 2), ...
    gmmpdfcalculator([xHat(nv1, 1) xHat(nv1, 2)], nuHat_g, sigmaHat_g, ...
    proportionHat_g), 'g.','markersize',12);
% hold on;
plot3(xHat(nv2, 1), xHat(nv2, 2), ...
    gmmpdfcalculator([xHat(nv2, 1) xHat(nv2, 2)], nuHat_g, sigmaHat_g, ...
    proportionHat_g), 'b.','markersize',12);

plot3(nuHat_k(:, 1), nuHat_k(:, 2), ...
    gmmpdfcalculator([nuHat_k(:, 1) nuHat_k(:, 2)], nuHat_g, sigmaHat_g, ...
    proportionHat_g), 'ro','markersize',12,'linewidth',8);

plot3(nuStar(:, 1), nuStar(:, 2), ...
    gmmpdfcalculator([nuStar(:, 1) nuStar(:, 2)], nuHat_g, sigmaHat_g, ...
    proportionHat_g), 'k+','markersize',18,'linewidth',8);

plot3(nuHat_g(:, 1), nuHat_g(:, 2), ...
    gmmpdfcalculator([nuHat_g(:, 1) nuHat_g(:, 2)], nuHat_g, sigmaHat_g, ...
    proportionHat_g), 'yx','markersize',18,'linewidth',8);

legend('1','2','k','truth','g')

%% Plot PDF
gridSize = 0.01;
[x, y] = meshgrid(-1:gridSize:1, -1:gridSize:1);

l = size(x, 1);
pointTmp = [reshape(x, l^2, 1), reshape(y, l^2, 1)];

% pdfGaussian = zeros(l^2, 1);
% for k = 1:2
%     pdfGaussian = pdfGaussian + proportionHat_g(k)*...
%         mvnpdf(pointTmp, nuHat_g(k, :), sigmaHat_g(:, :, k));
% end

pdfGaussian = gmmpdfcalculator(pointTmp, nuHat_g, sigmaHat_g, ...
    proportionHat_g);

z = reshape(pdfGaussian, l, l);
plot3(x, y, z, '-k');
view(45,45);








% nv1 = (tauStar == 1);
% nv2 = (tauStar == 2);
% 
% % axis([-1 1 -1 1 -1 1])
% hold on;
% plot(xHat(nv1, 1), xHat(nv1, 2), 'g.','markersize',12);
% hold on;
% plot(xHat(nv2, 1), xHat(nv2, 2), 'b.','markersize',12);
% 
% plot(nuHat_g(:, 1), nuHat_g(:, 2), 'ro','markersize',12,'linewidth',8);
% 
% if nargin>3
%     plot(nuStar(:, 1), nuStar(:, 2), 'k+','markersize',18,'linewidth',8);
% end
% 
% if nargin>4
%     plot(muHat2(:, 1), muHat2(:, 2), 'yx','markersize',18,'linewidth',8);
% end
% 
% legend('1','2','k','truth','g')