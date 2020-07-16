function [irf] = IRF_SVAR(By,ShockVector,nhorizons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(By,1);
nlags = size(By,3);
irf = zeros(nv,nhorizons + 1);
irf(:,1) = ShockVector;

for istep = 2:(nhorizons + 1)
    nlagsToUse = min(istep - 1,nlags);
    for ilag = 1:nlagsToUse
        irf(:,istep) = irf(:,istep) + By(:,:,ilag) * irf(:,istep - ilag);
    end
end

end

