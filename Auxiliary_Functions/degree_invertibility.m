function R0_2 = degree_invertibility(D, innov_var, shock_weight)

% Compute degree of invertibility of shock_weight'*e_t in state space model:
    %   s_t = A * s_{t-1} + B * e_t
    %   y_t = C * s_{t-1} + D * e_t
% with e_t ~ WN(0,I)

% Inputs:
% D             D matrix
% innov_var     var-cov matrix of Wold innovations u_t
% shock_weight  vector of shock weights with norm 1

aux = D*shock_weight; % cov(u_t, shock_weight'*e_t) = cov(y_t, shock_weight'*e_t)
R0_2 = aux'*(innov_var\aux); % degree of invertibility

end