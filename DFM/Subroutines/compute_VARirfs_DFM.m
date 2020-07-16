function IRF = compute_VARirfs_DFM(model, settings)

%   state-space form (A,B,C,D) based on Fernandez et al. (2005)
%   y_t^* = C * x_t + D w_t
%   x_{t+1} = A * x_t + B w_t
%   where y_t^* = (1 - delta(L)) y_t
%         x_t = (f_t, f_{t-1}, f_{t-2})
%   and w_t = (eta_{t+1}, v_t)

% prepare

IRF_hor = settings.est.IRF_hor;
n_fac = model.n_fac;
n_lags_fac = model.n_lags_fac;
n_lags_uar = model.n_lags_uar;

Phi = model.Phi;
Sigma_eta = model.Sigma_eta;
Lambda = model.Lambda;
delta = model.delta;
sigma_v = model.sigma_v;

var_select = settings.specifications.var_select;
n_spec = size(var_select,1);
n_var = size(var_select,2);
n_lags_state = max(n_lags_uar + 1, n_lags_fac);

recursive_shock_pos = settings.est.recursive_shock_pos;
IRF_response_var_pos = settings.est.IRF_response_var_pos;

IRF = NaN(IRF_hor, n_spec);

% compute A

A = zeros(n_lags_state * n_fac);
A(1:n_fac, 1:(n_lags_fac * n_fac)) = Phi(1:n_fac, :);
A((1 + n_fac):(n_lags_state * n_fac), 1:((n_lags_state - 1) * n_fac)) = eye((n_lags_state - 1) * n_fac);

% compute B

B = zeros(n_lags_state * n_fac, n_fac + n_var);
B(1:n_fac, 1:n_fac) = chol(Sigma_eta, 'lower');

% in each specification
for i_spec = 1:n_spec
    
    IRF_y_star = NaN(IRF_hor, n_var);
    IRF_y = NaN(IRF_hor, n_var);

    % compute C
    C_right = zeros(n_var * (n_lags_uar + 1), n_fac * n_lags_state);
    C_right(:, 1:(n_fac * (n_lags_uar + 1))) = kron(eye(n_lags_uar + 1), Lambda(var_select(i_spec, :), :));
    C_left = zeros(n_var, (n_lags_uar + 1) * n_var);
    C_left(:, 1:n_var) = eye(n_var);
    for ilag = 1:n_lags_uar
        C_left(:, ilag * n_var + (1:n_var)) = -diag(delta(var_select(i_spec, :), ilag));
    end
    C = C_left * C_right;
    
    % compute D

    D = zeros(n_var, n_fac + n_var);
    D(:, (n_fac + 1):end) = diag(sigma_v(var_select(i_spec, :), 1));
    
    % compute steady state conditional variance in Kalman filter
    cond_var = cond_var_fn_St_1(A, B, C, D);
    
    % compute cov-var matrix of innovations
    Sigma_innovation = C * cond_var * C' + D * D';
    
    % compute Kalman gain
    K = (A * cond_var * C' + B * D') / Sigma_innovation;
    
    % compute cholesky decomposition
    G = chol(Sigma_innovation, 'lower');
    
    % compute IR of y_star to innovation
    IRF_y_star(1, :) = G(:, recursive_shock_pos)';
    for i_hor = 1:(IRF_hor - 1)
        IRF_y_star(i_hor + 1, :) = (C * A^(i_hor - 1) * K * G(:, recursive_shock_pos))';
    end
    
    % compute IR of y to innovation
    for i_y = 1:n_var
        deltaFilter = delta(var_select(i_spec, :), :);
        IRF_y(:, i_y) = filter(1, [1, -deltaFilter(i_y, :)], IRF_y_star(:, i_y));
    end
    
    % store normalized IRF
    IRF(:, i_spec) = IRF_y(:, IRF_response_var_pos) / IRF_y(1, recursive_shock_pos);

end

end