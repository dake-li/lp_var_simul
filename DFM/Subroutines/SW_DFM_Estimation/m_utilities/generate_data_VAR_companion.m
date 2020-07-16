function [y_mat,shocks] = generate_data_VAR_companion(T,nrep,Q,M,G,n_initial)
% Model is 
% y(t) = Q*x(t)
% x(t) = M*x(t-1) + G*u(t)
%
% var(u) = eye
%
% Generate nrep replications of vectors of y of length T
%
% initial value of x is zero
% n_initial periods are generated to approximate stationary distribution 
%
nx = size(M,1);
nu = size(G,2);
ny = size(Q,1);
x = zeros(nx,nrep);
shocks = randn(nu,nrep,n_initial + T);
for t = 1:n_initial;
    x = M*x + G*shocks(:,:,t);
end;
y_mat = NaN(ny,nrep,T);
for t = 1:T;
    x = M*x + G*shocks(:,:,t + n_initial);
    y = Q*x;
    y_mat(:,:,t) = y;
end;
shocks = permute(shocks,[3,1,2]);
y_mat = permute(y_mat,[3,1,2]);
    
end

