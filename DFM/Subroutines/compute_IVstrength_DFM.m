function IV_strength = compute_IVstrength_DFM(model, settings)
% Function for computing population measure of IV strength for each DGP,
% i.e. the fraction of y_t being explained by z_t, after controlling lagged
% z_t and \bar{w}_t
    % first, incorporate the DGP of z_t into the DFM model to get an augmented DFM:
        % factor transition: f^*_t = \Phi^*(L) f^*_{t-1} + H \epsilon^*_t
        % observables:       w_t = \Lambda^* f^*_t + v^*_t
        % measurement error: v^*_{it} = \Delta^*_i(L) v^*_{i,t-1} + \Xi^*_i \xi^*_{it}
        
        % where w_t = (z_t, \bar{w}_t')'
        %       f^*_t = (f_t', z_t)'
        %       \epsilon^*_t = (\epsilon_t', \nu_t)
        % (See Technical Companion Note for details)
        
    % second, transform the augmented DFM into ABCD representation:
        % state transition:  s_{t+1} = A * s_t + B * \zeta_t
        % measurement eq:    w^*_t = C * s_t + D * \zeta_t

        % where w^*_t = (1 - \Delta(L)) * w_t
        % (See details in Technical Companion Note)
        
    % third, transform the ABCD representation into VAR(\infty) representation:
        % VAR(\infty): % w^*_t = \Psi^*(L) w^*_{t-1} + u_t (See details in Fernandez et al., 2005)
        % Use Var(u_t) to derive Var(y_t | lagged w_t)
    
    % fourth, compute the explained variance of y_t by IV:
        % Var(y_t | lagged w_t) - Var(y_t | z_t, lagged w_t) = IRF(y,h=0)^2 * \alpha^2 / (\alpha^2 + \sigma^2_\nu)
        
    % finally, compute the explained fraction of y_t by IV as:
        % 1 - Var(y_t | z_t, lagged w_t) / Var(y_t | lagged w_t)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

% unpack settings

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
% Rewrite Coefficient Matrices for Augmented DFM
%----------------------------------------------------------------

% original coefficient matrices in DFM

Phi = model.Phi;
Sigma_eta = model.Sigma_eta;
G = chol(Sigma_eta, 'lower');
Lambda = model.Lambda;
delta = model.delta;
sigma_v = model.sigma_v;
n_y = size(Lambda, 1);

% add z_t into augmented DFM to construct auxiliary form

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

% add z_t into each DGP

var_select = [var_select, (n_y + 1) * ones(n_spec, 1)];

% update dimension
n_fac = size(Lambda_aux, 2);
n_var = size(var_select, 2);

%----------------------------------------------------------------
% Represent Augmented DFM as ABCD
%----------------------------------------------------------------

% compute A

n_lags_state = max(n_lags_uar + 1, n_lags_fac);
A = zeros(n_lags_state * n_fac);
A(1:n_fac, 1:(n_lags_fac * n_fac)) = Phi_aux(1:n_fac, :);
A((1 + n_fac):(n_lags_state * n_fac), 1:((n_lags_state - 1) * n_fac)) = eye((n_lags_state - 1) * n_fac);

% compute B

B = zeros(n_lags_state * n_fac, n_fac + n_var);
B(1:n_fac, 1:n_fac) = G_aux;

% in each DGP
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
    % Represent ABCD as VAR(\infty)
    %----------------------------------------------------------------
    
    % compute steady state conditional variance in Kalman filter
    cond_var = cond_var_fn_St_1(A, B, C, D);
    
    % compute cov-var matrix of innovations
    Sigma_innovation = C * cond_var * C' + D * D';
    
    %----------------------------------------------------------------
    % IV-Explained and Total Variance of y_t with Lagged Controls
    %----------------------------------------------------------------
    
    % variance in projection-error
    
    variance_projection_error = Sigma_innovation(IV_est_normalize_var_pos, IV_est_normalize_var_pos);
    
    % variance explained by IV
    
    variance_explained_IV = IRF_initial^2 * IV_alpha^2 / (IV_alpha^2 + IV_sigma_v^2);
    
    % population IV strength
    
    IV_strength(i_spec) = variance_explained_IV / variance_projection_error;

end

end

