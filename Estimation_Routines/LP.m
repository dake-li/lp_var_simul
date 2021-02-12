function [Bc,Br,Bx,By,Sigma,Sxx] = LP(Y,recurShock,respV,nlags,horizon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)

% generate data matrices
[y,x,w] = LP_gen_data(Y,recurShock,respV,nlags,horizon);
X = [ones(size(x,1),1), x, w];
[Beta,Sigma,Sxx,~] = LS(y,X);

Bc = Beta(1); % constant
if recurShock == 1
    Br = []; % contemperaneous control
else
    Br = Beta(3:(recurShock+1));
end
Bx = Beta(2); % endogenous variable related to shock
By = Beta((recurShock + 2):end);
By = reshape(By,[nv,nlags]);

end

