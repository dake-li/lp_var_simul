function R0_sq = compute_invert_DFM(model,settings);
% Function for computing degree of invertibility of the true shock using
% the selected endogenous variables in one DGP
% (for observed-shock and IV experiments)

    % First, need to transform the selected DGP from DFM representation to a general ABCD
    % representation.
    
    % DFM representation:
        % factor transition: f_t = \Phi(L) f_{t-1} + H \epsilon_t
        % observables:       \bar{w}_t = \bar{\Lambda} f_t + v_t
        % measurement error: v_{it} = \Delta_i(L) v_{i,t-1} + \Xi_i \xi_{it}
    
    % ABCD representation:
        % state transition:  s_t = A * s_{t-1} + B * \epsilon_t
        % measurement eq:    \bar{w}^*_t = C * s_t + D * \xi_t
        
        % where \bar{w}^*_t = (I - \bar{\Delta}(L)) * \bar{w}_t
        %       s_t = (\epsilon_t', \epsilon_{t-1}', ...)'
        
    % Transforming formula can be found in Technical Companion Note
        
    % Then, degree of invertibility:
        % first compute Var(s_t | \bar{w}^*_t, \bar{w}^*_{t-1}, ...) using Kalman filter
        % then derive R0_sq = Var(shock_weight' * \epsilon_t | \bar{w}_t, \bar{w}_{t-1}, ...)

% unpack settings

n_fac = model.n_fac;
n_lags_fac = model.n_lags_fac;
n_lags_uar = model.n_lags_uar;

Phi = model.Phi;
Sigma_eta = model.Sigma_eta;
Lambda = model.Lambda;
delta = model.delta;
sigma_v = model.sigma_v;

var_select = settings.specifications.var_select;
n_spec = size(var_select, 1);
n_var = size(var_select, 2);

VMA_nlags = settings.est.VMA_nlags;
shock_weight = settings.est.shock_weight;

% compute (1 - \Phi(L))^{-1} polynomial

Phi_polynomial = zeros(n_fac, n_fac, VMA_nlags + 1); % (1 - \Phi(L))^(-1)
Phi_lag = eye(n_fac * n_lags_fac); % Big_Phi_Mat^lag
Phi_polynomial(:,:,1) = Phi_lag(1:n_fac, 1:n_fac);

G = chol(Sigma_eta, 'lower'); % Warning: corresponds to matrix H in our paper
Phi_polynomial_G = zeros(n_fac, n_fac, VMA_nlags + 1);
Phi_polynomial_G(:,:,1) = Phi_polynomial(:,:,1) * G;

for ilag = 1:VMA_nlags
    
    Phi_lag = Phi_lag * Phi;
    Phi_polynomial(:,:,ilag+1) = Phi_lag(1:n_fac, 1:n_fac);
    Phi_polynomial_G(:,:,ilag+1) = Phi_polynomial(:,:,ilag+1) * G;
    
end

Phi_polynomial_G_mat = reshape(Phi_polynomial_G,[n_fac, n_fac * (VMA_nlags+1)]); % (1 - \Phi(L))^{-1} * G

% placeholder for R0_sq

R0_sq = zeros(n_spec, 1);

% compute matrix A

A = zeros(n_fac * (1 + VMA_nlags));
A((n_fac+1):end, 1:(VMA_nlags*n_fac)) = eye(VMA_nlags * n_fac);
A = sparse(A);

% compute matrix B

B = [eye(n_fac); zeros(VMA_nlags * n_fac, n_fac)];
B = sparse(B);

% go thru each DGP

for i_spec = 1:n_spec
    
    % compute matrix C
    
    Lambda_select = Lambda(var_select(i_spec,:),:);
    delta_select = delta(var_select(i_spec,:),:);
    
    C0 = Lambda_select * Phi_polynomial_G_mat; % \Lambda * (1 - \Phi(L))^{-1} * G
    C = C0;
    
    % (1 - \Delta(L)) * \Lambda * (1 - \Phi(L))^{-1} * G
    
    for ilag = 1:n_lags_uar
        
        C_lag = -delta_select(:,ilag) .* C0;
        C_lag_truncate = C_lag(:, 1:((VMA_nlags + 1 - ilag) * n_fac));
        C_lag_shift = [zeros(n_var,n_fac * ilag), C_lag_truncate];
        C = C + C_lag_shift;
        
    end
    
    % compute matrix D
    
    sigma_v_select = sigma_v(var_select(i_spec,:),:);
    D = diag(sigma_v_select);
    
    % run Kalman filter to compute R0_sq
    
    cond_var_include_lag = cond_var_fn_St(A,B,C,D); % Var(s_t | \bar{w}^*_t, \bar{w}^*_{t-1}, ...)
    cond_var = cond_var_include_lag(1:n_fac, 1:n_fac); % Var(\epsilon_t | \bar{w}^*_t, \bar{w}^*_{t-1}, ...)
    R0_sq(i_spec, 1) = 1 - shock_weight' * cond_var * shock_weight; % Var(shock_weight' * \epsilon_t | \bar{w}^*_t, \bar{w}^*_{t-1}, ...)
    % equivalent to Var(shock_weight' * \epsilon_t | \bar{w}_t, \bar{w}_{t-1}, ...)
end

end