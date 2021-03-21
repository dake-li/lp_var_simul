function IRF = compute_VARirfs_DFM(model, settings)
% Function for computing true IRF for each DGP in recursive experiments

    % First, need to transform the selected DGP from DFM representation to a general ABCD
    % representation.
    
    % DFM representation:
        % factor transition: f_t = \Phi(L) f_{t-1} + H \epsilon_t
        % observables:       \bar{w}_t = \bar{\Lambda} f_t + v_t
        % measurement error: v_{it} = \Delta_i(L) v_{i,t-1} + \Xi_i \xi_{it}
    
    % ABCD representation:
        % state transition:  s_{t+1} = A * s_t + B * \zeta_t
        % measurement eq:    w^*_t = C * s_t + D * \zeta_t
        
        % where w^*_t = (I - \Delta(L)) * w_t
        %       s_t = (f_t', f_{t-1}', f_{t-2}')'
        %       \zeta_t = (\epsilon_t', \xi_t')'
        
    % Transforming formula can be found in Technical Companion Note
        
    % Then, transform this ABCD representation into VAR(\infty) representation:
        % w^*_t = \Psi^*(L) w^*_{t-1} + u_t (See details in Fernandez et al., 2005)
    
    % Finally, in VAR(\infty) representation:
        % first derive IRF of w^*_t to the recursive shock
        % and then derive IRF of w_t to the recursive shock

%------------------------------
% Prepare
%------------------------------
        
% unpack settings

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

%------------------------------
% Transform DFM to ABCD
%------------------------------

% compute matrix A

A = zeros(n_lags_state * n_fac);
A(1:n_fac, 1:(n_lags_fac * n_fac)) = Phi(1:n_fac, :);
A((1 + n_fac):(n_lags_state * n_fac), 1:((n_lags_state - 1) * n_fac)) = eye((n_lags_state - 1) * n_fac);

% compute matrix B

B = zeros(n_lags_state * n_fac, n_fac + n_var);
B(1:n_fac, 1:n_fac) = chol(Sigma_eta, 'lower');

% in each DGP
for i_spec = 1:n_spec
    
    IRF_y_star = NaN(IRF_hor, n_var); % Warning: y_star corresponds to w^* in our paper
    IRF_y = NaN(IRF_hor, n_var);

    % compute matrix C
    C_right = zeros(n_var * (n_lags_uar + 1), n_fac * n_lags_state);
    C_right(:, 1:(n_fac * (n_lags_uar + 1))) = kron(eye(n_lags_uar + 1), Lambda(var_select(i_spec, :), :));
    C_left = zeros(n_var, (n_lags_uar + 1) * n_var);
    C_left(:, 1:n_var) = eye(n_var);
    for ilag = 1:n_lags_uar
        C_left(:, ilag * n_var + (1:n_var)) = -diag(delta(var_select(i_spec, :), ilag));
    end
    C = C_left * C_right;
    
    % compute matrix D

    D = zeros(n_var, n_fac + n_var);
    D(:, (n_fac + 1):end) = diag(sigma_v(var_select(i_spec, :), 1));
    
    %------------------------------
    % Transform ABCD to VAR(\infty)
    %------------------------------
    
    % compute steady state conditional variance in Kalman filter
    cond_var = cond_var_fn_St_1(A, B, C, D);
    
    % compute cov-var matrix of innovations
    Sigma_innovation = C * cond_var * C' + D * D';
    
    % compute Kalman gain
    K = (A * cond_var * C' + B * D') / Sigma_innovation;
    
    % compute cholesky decomposition
    G = chol(Sigma_innovation, 'lower');
    
    %------------------------------
    % Compute IRF in VAR(\infty)
    %------------------------------
    
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