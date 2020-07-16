function [Bc,Br,Bx,By,Sigma,Sxx] = LP(Y,recurShock,respV,nlags,horizon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)
    
X = lagmatrix(Y,horizon + (1:nlags));
X = [lagmatrix(Y(:, 1:recurShock), horizon), X];
Y = Y((horizon + nlags + 1):end, respV);
X = X((horizon + nlags + 1):end, :);
X = [ones(size(X,1),1), X];
[Beta,Sigma,Sxx,~] = LS(Y,X);

Bc = Beta(1,1);
if recurShock == 1
    Br = [];
else
    Br = Beta(2:recurShock, 1);
end
Bx = Beta(recurShock + 1, 1);
By = Beta((recurShock + 2):end, 1);
By = reshape(By,[nv,nlags]);

end

