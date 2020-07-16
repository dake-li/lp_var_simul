function irf = IRF_LP(Y,recurShock,respV,nlags,nhorizons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2);
nT = size(Y,1);
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end
irf = zeros(1,nhorizons + 1);

for h = 0:nhorizons
    [~,~,Bx,~,~,~] = LP(Y,recurShock,respV,nlags,h);
    irf(1,h + 1) = Bx;
end

end

