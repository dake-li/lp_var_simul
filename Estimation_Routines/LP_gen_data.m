function [y,x,w] = LP_gen_data(Y,recurShock,respV,nlags,nhorizon)
% Function for generating standard data matrices for LP

nv = size(Y,2);
nT = size(Y,1);

y  = Y(:,respV); % response variable
x  = lagmatrix(Y(:,recurShock), nhorizon); % impulse variable
if recurShock > 1 % check if there are contemperaneous controls
    w  = [ lagmatrix(Y(:,1:(recurShock - 1)), nhorizon) , lagmatrix( Y , (1:nlags) + nhorizon ) ]; % control variables (contemporaneous vars, lagged vars, no constant)
else
    w  = lagmatrix( Y , (1:nlags) + nhorizon );
end
y = y((nlags + nhorizon + 1):end);
x = x((nlags + nhorizon + 1):end);
w = w((nlags + nhorizon + 1):end, :); % Warning: here w include both contemperaneous and lagged controls
end

