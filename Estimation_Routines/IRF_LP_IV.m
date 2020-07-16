function [irf, Fstat_z] = IRF_LP_IV(H,respV,normlzV,nlags,nhorizons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(H,2) - 1;
nT = size(H,1);
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end
irf = zeros(1, nhorizons + 1);

for h = 0:nhorizons
    [firstStage,secondStage] = LP_IV(H,respV,normlzV,nlags,h);
    irf(1,h + 1) = secondStage.Bz / firstStage.Bz;
end

Fstat_z = firstStage.Fstat_z;

end

