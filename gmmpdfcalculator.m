function pdfGaussian = gmmpdfcalculator(pointTmp, nuHat_g, sigmaHat_g, ...
    proportionHat_g)

pdfGaussian = zeros(size(pointTmp, 1), 1);
for k = 1:size(sigmaHat_g, 3)
    pdfGaussian = pdfGaussian + proportionHat_g(k)*...
        mvnpdf(pointTmp, nuHat_g(k, :), sigmaHat_g(:, :, k));
end
