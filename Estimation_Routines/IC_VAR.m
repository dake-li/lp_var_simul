function [BIC,AIC] = IC_VAR(Y,nlagsMax)
% Function for computing BIC and AIC in choosing VAR lag length

nv = size(Y,2);
nT = size(Y,1);
BIC = zeros(1,nlagsMax);
AIC = zeros(1,nlagsMax);

% go thru from one lag to max number of lags
for nlags = 1:nlagsMax
    [~,~,Sigma,~,~] = VAR(Y((nlagsMax - nlags + 1):end,:),nlags);
    BIC(1,nlags) = log(det(Sigma)) + (nv^2 * nlags + nv) * log(nT - nlagsMax) / (nT - nlagsMax);
    AIC(1,nlags) = log(det(Sigma)) + (nv^2 * nlags + nv) * 2 / (nT - nlagsMax);
end

end

