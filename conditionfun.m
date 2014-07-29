function [c, ceq] = conditionfun(x, nVertex, dimLatentPosition)
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
xHat = reshape(x, nVertex, dimLatentPosition);

%% Inequalities
c = [];
for iVertex = 1:nVertex
    c = [c; xHat(iVertex, :)*xHat(iVertex, :)' - 1];
end

%% Equalities
ceq = 0;

end