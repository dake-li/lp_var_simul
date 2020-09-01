function cond_var = cond_var_fn_St_1(A, B, C, D)
%UNTITLED Summary of this function goes here
%   state-space form (A,B,C,D) based on FRS paper (2005)
%   y_t = C * x_t + D w_t
%   x_{t+1} = A * x_t + B w_t
%   compute var-cov of x_t in innovation representation

% number of state variables
n_s = size(C,2);

% initial conditions
Pt = eye(n_s);

% stop condition
dist = 1;
tol = 10^(-10);
relax = 0.9;

% Kalman filtering
while dist >= tol
    
    % prediction equations
    Pt_1 = A * Pt * A' + B * B';
    
    % updating equations
    ht = C * Pt * C' + D * D';
    qt = A * Pt * C' + B * D';
    Pt_upd = Pt_1 - qt / ht * qt';
    
    % update P_t
    dist = max(max(abs(Pt_upd - Pt)));
    Pt = relax * Pt + (1 - relax) * Pt_upd;
    
end

cond_var = Pt;

end

