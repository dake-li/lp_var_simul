function [y,x,w] = LP_gen_data(Y,recurShock,respV,nlags,nhorizon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)
nT = size(Y,1);

y  = Y(:,respV); % endogenous variable
x  = lagmatrix(Y(:,recurShock), nhorizon); % endoegnous variable related to the shock
if recurShock > 1
    w  = [ lagmatrix(Y(:,1:(recurShock - 1)), nhorizon) , lagmatrix( Y , (1:nlags) + nhorizon ) ]; % control variables (contemporaneous vars, lagged vars, no constant)
else
    w  = lagmatrix( Y , (1:nlags) + nhorizon );
end
y = y((nlags + nhorizon + 1):end);
x = x((nlags + nhorizon + 1):end);
w = w((nlags + nhorizon + 1):end, :);
end

