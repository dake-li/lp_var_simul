function [BIC,AIC] = IC_LP(Y,dimR,nlagsMax,horizon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2);
nT = size(Y,1);
BIC = zeros(1,nlagsMax);
AIC = zeros(1,nlagsMax);

for nlags = 1:nlagsMax
    [Bc,Br,Bx,By,Sigma,Sxx] = LP(Y((nlagsMax - nlags + 1):end,:),dimR,nlags,horizon);
    BIC(1,nlags) = log(det(Sigma)) + (dimR + 1 + nv * nlags) * log(nT - horizon - nlagsMax) / (nT - horizon - nlagsMax);
    AIC(1,nlags) = log(det(Sigma)) + (dimR + 1 + nv * nlags) * 2 / (nT - horizon - nlagsMax);
end

end

