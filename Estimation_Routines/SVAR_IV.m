function [VARout,IVout] = SVAR_IV(H,normlzV,nlags)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% put the shock of interest at the same position as x
Z = H(:,1); % h = (z,r,x,y,q)
Y = H(:,2:end);
nv = size(Y,2);
nT = size(Y,1);

% IV
X = lagmatrix(Y,1:nlags);
X = [ones(nT,1),Z,X];
Y = Y((nlags + 1):end,:);
X = X((nlags + 1):end,:);

[Beta,Sigma,Sxx,~] = LS(Y,X);
Theta = Beta(2,:)';
IVout.gamma = Theta / Theta(normlzV); % normalize
IVout.Sigma = Sigma;
IVout.Sxx = Sxx;
Avar = Sigma(normlzV, normlzV) * inv(Sxx);
IVout.Fstat_z = size(Y, 1) * Theta(normlzV)^2 / Avar(2, 2); % Wald stat

% VAR
Y = H(:,2:end);
[VARout.Bc,VARout.By,VARout.Sigma,VARout.Sxx] = VAR(Y,nlags);

end