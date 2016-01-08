

nVec = [150, 300, 600, 900];
nN = length(nVec);


timeEBSBM = zeros(1, nN);

for iN = 1:nN
    n = nVec(iN);
    timeEBSBM(iN) = ebsbmsim(n, 3, 0.1, 1, 50, 2, 2, 2, 2);
end


timeOpti = zeros(1, nN);

for iN = 1:nN
    n = nVec(iN);
    timeOpti(iN) = sbmoptisim(n, 3, 0.3, 0.1, -1, 1, 50, 0);
end