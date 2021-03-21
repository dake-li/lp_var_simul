function cond_var = cond_var_fn_St_1(A, B, C, D)
% Function to compute the conditional variance, Var(s_t | y_t, y_{t-1},...) using the Kalman filter
    % This function applies to the state-space form of the following form:
        %   s_{t+1} = A * s_t + B w_t
        %   y_t = C * s_t + D w_t
        % (state-space form based on Fernandez et al., 2005)

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

