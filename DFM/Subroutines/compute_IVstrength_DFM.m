function IV_strength = compute_IVstrength_DFM(model, settings)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_fac = model.n_fac;
n_lags_fac = model.n_lags_fac;
n_lags_uar = model.n_lags_uar;

IV_rho = model.IV.rho;
IV_alpha = model.IV.alpha;
IV_sigma_v = model.IV.sigma_v;

var_select = settings.specifications.var_select;
n_spec = size(var_select,1);
n_var = size(var_select,2);

shock_weight = settings.est.shock_weight;
IV_est_normalize_var_pos = settings.est.IV_est_normalize_var_pos;

IV_strength = NaN(n_spec, 1);

%----------------------------------------------------------------
% Rewrite Coefficient Matrices for DFM
%----------------------------------------------------------------

% original coefficient matrices in DFM

Phi = model.Phi;
Sigma_eta = model.Sigma_eta;
G = chol(Sigma_eta, 'lower');
Lambda = model.Lambda;
delta = model.delta;
sigma_v = model.sigma_v;
n_y = size(Lambda, 1);

% add z_t into DFM to construct auxiliary form

Phi_aux = zeros((n_fac + 1) * n_lags_fac);
for i_lag_fac = 1:n_lags_fac
    Phi_aux(1:n_fac, (i_lag_fac - 1) * (n_fac + 1) + (1:n_fac)) = Phi(1:n_fac, (i_lag_fac - 1) * n_fac + (1:n_fac));
end
Phi_aux(n_fac + 1, n_fac + 1) = IV_rho;
if n_lags_fac > 1
    Phi_aux((n_fac + 2):end, 1:((n_lags_fac - 1) * (n_fac + 1))) = eye((n_lags_fac - 1) * (n_fac + 1));
end

G_aux = [G, zeros(n_fac, 1); zeros(1, n_fac), IV_sigma_v];
G_aux(n_fac + 1, 1:n_fac) = IV_alpha * shock_weight';

Lambda_aux = [Lambda, zeros(n_y, 1); zeros(1, n_fac), 1];

delta_aux = [delta; zeros(1, n_lags_uar)];
sigma_v_aux = [sigma_v; 0];

% add z_t into each specification

var_select = [var_select, (n_y + 1) * ones(n_spec, 1)];

% update dimension
n_fac = size(Lambda_aux, 2);
n_var = size(var_select, 2);

%----------------------------------------------------------------
% represent DFM as ABCD
%----------------------------------------------------------------

% compute A

n_lags_state = max(n_lags_uar + 1, n_lags_fac);
A = zeros(n_lags_state * n_fac);
A(1:n_fac, 1:(n_lags_fac * n_fac)) = Phi_aux(1:n_fac, :);
A((1 + n_fac):(n_lags_state * n_fac), 1:((n_lags_state - 1) * n_fac)) = eye((n_lags_state - 1) * n_fac);

% compute B

B = zeros(n_lags_state * n_fac, n_fac + n_var);
B(1:n_fac, 1:n_fac) = G_aux;

% in each specification
for i_spec = 1:n_spec
    
    % IRF of normalization variable at h = 0
    
    IRF_initial = model.irf(1,var_select(i_spec, IV_est_normalize_var_pos));

    % compute C
    C_right = zeros(n_var * (n_lags_uar + 1), n_fac * n_lags_state);
    C_right(:, 1:(n_fac * (n_lags_uar + 1))) = kron(eye(n_lags_uar + 1), Lambda_aux(var_select(i_spec, :), :));
    C_left = zeros(n_var, (n_lags_uar + 1) * n_var);
    C_left(:, 1:n_var) = eye(n_var);
    for ilag = 1:n_lags_uar
        C_left(:, ilag * n_var + (1:n_var)) = -diag(delta_aux(var_select(i_spec, :), ilag));
    end
    C = C_left * C_right;
    
    % compute D

    D = zeros(n_var, n_fac + n_var);
    D(:, (n_fac + 1):end) = diag(sigma_v_aux(var_select(i_spec, :), 1));
    
    %----------------------------------------------------------------
    % Compute Variance of Reduced-form Errors
    %----------------------------------------------------------------
    
    % compute steady state conditional variance in Kalman filter
    cond_var = cond_var_fn_St_1(A, B, C, D);
    
    % compute cov-var matrix of innovations
    Sigma_innovation = C * cond_var * C' + D * D';
    
    % variance in projection-error
    
    variance_projection_error = Sigma_innovation(IV_est_normalize_var_pos, IV_est_normalize_var_pos);
    
    % variance explained by IV
    
    variance_explained_IV = IRF_initial^2 * IV_alpha^2 / (IV_alpha^2 + IV_sigma_v^2);
    
    % population IV strength
    
    IV_strength(i_spec) = variance_explained_IV / variance_projection_error;

end

end

