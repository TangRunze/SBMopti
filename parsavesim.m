function  parsavesim(fName, errorRateASGE, errorRateOpti, xBest, nuBest,...
    tauBest, xHat0, nuHat0, tauHat0, loglikASGE, loglikOpti);
% save data
save(fName, 'errorRateASGE', 'errorRateOpti', 'xBest', 'nuBest', ...
    'tauBest', 'xHat0', 'nuHat0', 'tauHat0', 'loglikASGE', 'loglikOpti');
end