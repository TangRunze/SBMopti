function val = projectobjectivefun(x, dimLatentPosition, xHatTmp)
% Objective functions in optimization problem

rotationMatrix = reshape(x, dimLatentPosition, dimLatentPosition);

xHat = xHatTmp*rotationMatrix;

val = - min(min(abs(xHat)));

end