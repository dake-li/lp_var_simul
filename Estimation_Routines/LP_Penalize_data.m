function [y,x,w,H_min,H_max,r] = LP_Penalize_data(Y,recurShock,respV,nlags,nhorizon,irfLimitOrder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)
nT = size(Y,1);

H_min = 0;
H_max = nhorizon;

y  = Y(:,respV); % endogenous variable
x  = Y(:,recurShock); % endoegnous variable related to the shock
if recurShock > 1
    w  = [ Y(:,1:(recurShock - 1)) , lagmatrix( Y , 1:nlags ) ]; % control variables (contemporaneous vars, lagged vars)
else
    w  = lagmatrix( Y , 1:nlags );
end
w( ~isfinite(w) ) = 0;

r = irfLimitOrder + 1;

end

