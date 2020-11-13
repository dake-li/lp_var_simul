function [Bc,By,Sigma,Sxx,Res,Beta,Y,X] = VAR(Y,nlags)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)
X = lagmatrix(Y,1:nlags);
Y = Y((nlags + 1):end,:);
X = X((nlags + 1):end,:);
X = [ones(size(X,1),1),X];
[Beta,Sigma,Sxx,Res] = LS(Y,X);

Bc = Beta(1,:)';
By = reshape(Beta(2:end,:),[nv,nlags,nv]);
By = permute(By,[3,1,2]); % nv * nv * nlags

end

