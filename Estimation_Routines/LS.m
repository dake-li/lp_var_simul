function [Beta,Sigma,Sxx,Res] = LS(Y,X)
% Function for OLS regression

% error checking
if size(Y,1)~=size(X,1)
    error("Number of obs in Y and X do not match!")
end
if size(X,1)<size(X,2)
    error("Too few obs in regression!")
end

% OLS estimate
Beta = X \ Y;
Res = Y - X * Beta;
Sigma = Res' * Res / (size(Y,1) - size(X, 2));
Sxx = X' * X / size(Y,1);

end

