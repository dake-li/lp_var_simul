function cond_var = cond_var_state(A, B, C, D)
% Function to compute the conditional variance Var(s_t | y_t, y_{t-1},...)
% for state space models of the following form:
    %   s_t = A * s_{t-1} + B * e_t
    %   y_t = C * s_{t-1} + D * e_t
% with e_t ~ WN(0,I)

% number of state variables and shocks
[n_s,n_e] = size(B);

% initial conditions
P = 100*eye(n_s);

% stop condition
tol = 1e-7; % convergence tolerance
relax = 0; % weight to put on past values to slow down updates

% Iteratively solve Riccati equation
% Eq. 9 in Fernandez-Villaverde, Rubio-Ramirez, Sargent & Watson (AER 2007)
% Here implemented so as to ensure positive semidefiniteness

dist = 1+tol;

while dist >= tol
    
    PI = blkdiag(P,eye(n_e));
    AB = [A B];
    CD = [C D];
    L = AB-AB*PI*CD'*((CD*PI*CD')\CD);

    P_upd = L*PI*L';
    
    % update
    dist = max((abs(P_upd(:) - P(:))));
    P = relax * P + (1 - relax) * P_upd;
    
end

cond_var = P;

end

