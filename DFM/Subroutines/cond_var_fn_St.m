function cond_var = cond_var_fn_St(A,B,C,D)
% Function to compute the conditional variance, Var(s_t | y_t, y_{t-1},...) using Kalman filter
    % This function applies to the state-space form of the following form:
        %   y(t) = C * s(t) + D * e(t)
        %   s(t) = A * s(t-1) + B * u(t)

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
    ht = C * Pt_1 * C' + D * D';
    Kt = Pt_1 * C' * ht^(-1);
    Pt_upd = Pt_1 - Kt * C * Pt_1;
    
    % update P_t
    dist = max(max(abs(Pt_upd - Pt)));
    Pt = relax * Pt + (1 - relax) * Pt_upd;
    
end

cond_var   = Pt;

end

