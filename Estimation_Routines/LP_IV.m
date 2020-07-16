function [firstStage,secondStage] = LP_IV(H,respV,normlzV,nlags,horizon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Z = H(:,1); % h = (z,r,x,y,q)
Y = H(:,2:end);
nv = size(Y,2);
nT = size(Y,1);

% first stage
X = lagmatrix(H,1:nlags);
X = [ones(nT,1),Z,X];
Y = Y((nlags + 1):end, normlzV);
X = X((nlags + 1):end,:);
[Beta,Sigma,Sxx,~] = LS(Y,X);

firstStage.Bc = Beta(1,1);
firstStage.Bz = Beta(2,1);
By = Beta(3:end,1);
firstStage.By = reshape(By,[nv+1,nlags]);
firstStage.Sigma = Sigma;
firstStage.Sxx = Sxx;
Avar = Sigma * inv(Sxx);
firstStage.Fstat_z = size(Y, 1) * Beta(2,1)^2 / Avar(2, 2);

% second stage
Y = H(:,2:end);
X = lagmatrix(H,horizon + (1:nlags));
X = [ones(nT,1),lagmatrix(Z,horizon),X];
Y = Y((horizon + nlags + 1):end, respV);
X = X((horizon + nlags + 1):end,:);
[Beta,Sigma,Sxx,~] = LS(Y,X);

secondStage.Bc = Beta(1,1);
secondStage.Bz = Beta(2,1);
By = Beta(3:end,1);
secondStage.By = reshape(By,[nv+1,nlags]);
secondStage.Sigma = Sigma;
secondStage.Sxx = Sxx;
    
end

