function [VARout,IVout] = SVAR_IV(H,normlzV,nlags)
% Auxiliary function for identifying a SVAR using an external IV series

% split IV and endogenous variables in the data matrix
Z = H(:,1);
Y = H(:,2:end); % Warning: H corresponds to w_t, Y corresponds to \bar{w}_t in our paper 
nv = size(Y,2);
nT = size(Y,1);

% put IV and lagged endogeneous variables into the regressors X
X = lagmatrix(Y,1:nlags);
X = [ones(nT,1),Z,X];
Y = Y((nlags + 1):end,:);
X = X((nlags + 1):end,:);

% run multiple regression for each endogenous variable for impact response (equivalent to multiple 2SLS with IV for normalization variable)
[Beta,Sigma,Sxx,~] = LS(Y,X);
Theta = Beta(2,:)'; % impact response identified by IV
IVout.gamma = Theta / Theta(normlzV); % normalized response
IVout.Sigma = Sigma;
IVout.Sxx = Sxx;
Avar = Sigma(normlzV, normlzV) * inv(Sxx);
IVout.Fstat_z = size(Y, 1) * Theta(normlzV)^2 / Avar(2, 2); % Wald stat to test IV strength

% run least-squares VAR to estimate VAR coeffcients for lagged terms
Y = H(:,2:end);
[VARout.Bc,VARout.By,VARout.Sigma,VARout.Sxx] = VAR(Y,nlags);

end