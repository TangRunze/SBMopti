function errorRate = errorratecalculator(tauStar, tauHat, nVertex, nBlock)
% Calculate error rate
errorRate = nVertex;
if (tauStar ~= 0)
    permutation = perms(1:nBlock);
    for iFactorial = 1:factorial(nBlock)
        position = permutation(iFactorial, :);
        tmpTau = tauHat;
        for jBlock = 1:nBlock
            nv = (tauHat == position(jBlock));
            tmpTau(nv) = jBlock;
        end
        if sum(tauStar ~= tmpTau) < errorRate
            errorRate = sum(tauStar ~= tmpTau);
        end
    end
end
errorRate = errorRate/nVertex;
end
