function [c, ceq] = projectconditionfun(x, nVertex, dimLatentPosition, ...
    xHatTmp)
% Condition functions in optimization problem

% %% Pre-calculation
% xHat = reshape(x, nVertex, dimLatentPosition);
% P = xHat*xHat';
% 
% %% Inequalities
% c = [];
% for iVertex = 1:nVertex
%     for jVertex = 1:nVertex
%         c = [c; - P(iVertex, jVertex); P(iVertex, jVertex) - 1];
%     end
% end

%% Pre-calculation
rotationMatrix = reshape(x, dimLatentPosition, dimLatentPosition);
xHatTmp = xHatTmp*rotationMatrix;

%% Inequalities
c = [];
for iVertex = 1:nVertex
    for jDim = 1:dimLatentPosition
        c = [c; - xHatTmp(iVertex, jDim)];
    end
end

%% Equalities
ceq = [];
for iDim = 1:dimLatentPosition
    for jDim = iDim:dimLatentPosition
        if (iDim == jDim)
            ceq = [ceq; rotationMatrix(iDim, :)*...
                rotationMatrix(jDim, :)' - 1];
        else
            ceq = [ceq; rotationMatrix(iDim, :)*...
                rotationMatrix(jDim, :)'];
        end
    end
end

end