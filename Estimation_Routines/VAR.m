function [Bc,By,Sigma,Sxx,Res,Beta,Y,X] = VAR(Y,nlags)
% Auxiliary function for estimating least-squares VAR coefficients

% prepare data matrices
nv = size(Y,2);
X = lagmatrix(Y,1:nlags); % regressors X
Y = Y((nlags + 1):end,:);
X = X((nlags + 1):end,:);
X = [ones(size(X,1),1),X];

% OLS
[Beta,Sigma,Sxx,Res] = LS(Y,X);

% store coefficients
Bc = Beta(1,:)'; % constant term
By = reshape(Beta(2:end,:),[nv,nlags,nv]); % coeffcients on lagged terms
By = permute(By,[3,1,2]); % dim: nv * nv * nlags

end

