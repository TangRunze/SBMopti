function loglik = loglikcalculator(adjMatrix, tauHat, nBlock)

loglik = 0;
tmpB = zeros(nBlock, nBlock);
nv = cell(1, nBlock);
for i = 1:nBlock
    nv{i} = (tauHat == i);
end
for i = 1:nBlock
    for j = 1:nBlock
        tmpSumB = sum(sum(adjMatrix(nv{i}, nv{j})));
        if (i == j)
            tmpB(i, j) = tmpSumB/sum(nv{i})/(sum(nv{i}) - 1);
            if (tmpB(i, j) > 0) && (tmpB(i, j) < 1)
                loglik = loglik + tmpSumB*log(tmpB(i, j)) + ...
                    (sum(nv{i})*(sum(nv{i})-1)-tmpSumB)*log(1-tmpB(i, j));
            end
        else
            tmpB(i, j) = tmpSumB/sum(nv{i})/sum(nv{j});
            if (tmpB(i, j) > 0) && (tmpB(i, j) < 1)
                loglik = loglik + tmpSumB*log(tmpB(i, j)) + ...
                    (sum(nv{i})*sum(nv{j})-tmpSumB)*log(1-tmpB(i, j));
            end
        end
    end
end

% loglik = 0;
% for i = 1:(nVertex - 1)
%     for j = (i+1):nVertex
%         tmpP = tmpB(tauHat(i), tauHat(j));
%         if (tmpP > 0)
%             loglik = loglik + adjMatrix(i, j)*log(tmpP) + ...
%                 (1 - adjMatrix(i, j))*log(1 - tmpP);
%         end
%     end
% end
% loglik = loglik*2;

