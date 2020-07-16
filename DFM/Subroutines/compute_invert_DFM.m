function R0_sq = compute_invert_DFM(model,settings);

%   state-space form (A,B,C,D)
%   y_t^* = C * eta_t^* + D v_t 
%         = (1 - delta(L)) Lambda (1 - Phi(L))^{-1} G eta_t^* + diag(sigma_v) v_t
%   eta_t^* = A * eta_{t-1}^* + B u_t 
%           = [0,0;1,0] eta_{t-1}^* + [1;0] u_t
%   where y_t^* = (1 - delta(L)) y_t
%   and eta_^* is dg_invert_n_lags lagged terms of eta_t

% preparation

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

% compute (1 - Phi(L))^{-1} polynomial

Phi_polynomial = zeros(n_fac, n_fac, VMA_nlags + 1); % (1 - Phi(L))^(-1)
Phi_lag = eye(n_fac * n_lags_fac); % Big_Phi_Mat^lag
Phi_polynomial(:,:,1) = Phi_lag(1:n_fac, 1:n_fac);

G = chol(Sigma_eta, 'lower');
Phi_polynomial_G = zeros(n_fac, n_fac, VMA_nlags + 1);
Phi_polynomial_G(:,:,1) = Phi_polynomial(:,:,1) * G;

for ilag = 1:VMA_nlags
    
    Phi_lag = Phi_lag * Phi;
    Phi_polynomial(:,:,ilag+1) = Phi_lag(1:n_fac, 1:n_fac);
    Phi_polynomial_G(:,:,ilag+1) = Phi_polynomial(:,:,ilag+1) * G;
    
end

Phi_polynomial_G_mat = reshape(Phi_polynomial_G,[n_fac, n_fac * (VMA_nlags+1)]);

% in each specification

R0_sq = zeros(n_spec, 1);

% compute A

A = zeros(n_fac * (1 + VMA_nlags));
A((n_fac+1):end, 1:(VMA_nlags*n_fac)) = eye(VMA_nlags * n_fac);
A = sparse(A);

% compute B

B = [eye(n_fac); zeros(VMA_nlags * n_fac, n_fac)];
B = sparse(B);

for i_spec = 1:n_spec
    
    % compute C
    
    Lambda_select = Lambda(var_select(i_spec,:),:);
    delta_select = delta(var_select(i_spec,:),:);
    
    C0 = Lambda_select * Phi_polynomial_G_mat;
    C = C0;
    
    for ilag = 1:n_lags_uar
        
        C_lag = -delta_select(:,ilag) .* C0;
        C_lag_truncate = C_lag(:, 1:((VMA_nlags + 1 - ilag) * n_fac));
        C_lag_shift = [zeros(n_var,n_fac * ilag), C_lag_truncate];
        C = C + C_lag_shift;
        
    end
    
    % compute D
    
    sigma_v_select = sigma_v(var_select(i_spec,:),:);
    D = diag(sigma_v_select);
    
    % run Kalman filter to compute R0_sq
    
    cond_var_include_lag = cond_var_fn_St(A,B,C,D);
    cond_var = cond_var_include_lag(1:n_fac, 1:n_fac);
    R0_sq(i_spec, 1) = 1 - shock_weight' * cond_var * shock_weight;
    
end

end