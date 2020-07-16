function IV_strength = compute_IVstrength_DSGE(model, settings)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

var_select = settings.specifications.var_select;
n_spec = size(var_select, 1);

shock_weight = settings.est.shock_weight;
IV_est_normalize_var_pos = settings.est.IV_est_normalize_var_pos;

IV_rho = model.IV.rho;
IV_alpha = model.IV.alpha;
IV_sigma_v = model.IV.sigma_v;

IV_strength = NaN(n_spec, 1);

%----------------------------------------------------------------
% Rewrite Coefficient Matrices for ABCD
%----------------------------------------------------------------

% original ABCD matrices

A   = model.ABCD.A;
B   = model.ABCD.B;
C   = model.ABCD.C;
D   = model.ABCD.D;
n_s = size(A, 1);
n_y = size(C, 1);
n_eps = size(B, 2);

% add z_t into ABCD to construct auxiliary form
    
A_aux = [A, zeros(n_s, 1); zeros(1, n_s), IV_rho];
B_aux = [B, zeros(n_s, 1); IV_alpha * shock_weight', IV_sigma_v];
C_aux = [C, zeros(n_y, 1); zeros(1, n_s), IV_rho];
D_aux = [D, zeros(n_y, 1); IV_alpha * shock_weight', IV_sigma_v];

% add z_t into each specification

var_select = [var_select, (n_y + 1) * ones(n_spec, 1)];

for i_spec = 1:n_spec
    
    % extract coefficient matrices
    
    A_mat   = A_aux;
    B_mat   = B_aux;
    C_mat   = C_aux(var_select(i_spec,:),:);
    D_mat   = D_aux(var_select(i_spec,:),:);
    n_s = size(A_mat, 1);
    n_y = size(C_mat, 1);
    
    % IRF of normalization variable at h = 0
    
    IRF_initial = model.irf(1,var_select(i_spec, IV_est_normalize_var_pos));
    
    %----------------------------------------------------------------
    % Innovation Representation
    %----------------------------------------------------------------

    % derive Sigma and K

    Sigma_states = eye(n_s);
    K_states = zeros(n_s,n_y);

    dist = 1;
    tol = 10^(-10);
    relax = 0.9;

    while dist >= tol
        Sigma_upd = (A_mat - K_states*C_mat)*Sigma_states*(A_mat-K_states*C_mat)' +...
            B_mat*B_mat' + K_states*(D_mat*D_mat')*K_states' -...
            B_mat * D_mat'*K_states' - K_states * D_mat * B_mat';
        K_upd = (A_mat * Sigma_upd * C_mat' + B_mat * D_mat') * (C_mat * Sigma_upd * C_mat' + D_mat * D_mat')^(-1);

        dist1 = max(max(abs(Sigma_upd - Sigma_states)));
        dist2 = max(max(abs(K_upd - K_states)));

        Sigma_states = relax * Sigma_states + (1-relax) * Sigma_upd;
        K_states = relax * K_states + (1-relax) * K_upd;

        dist = max(dist1, dist2);
    end
    
    % variance in projection-error
    
    Sigma_innovation = C_mat * Sigma_states * C_mat' + D_mat * D_mat';
    variance_projection_error = Sigma_innovation(IV_est_normalize_var_pos, IV_est_normalize_var_pos);
    
    % variance explained by IV
    
    variance_explained_IV = IRF_initial^2 * IV_alpha^2 / (IV_alpha^2 + IV_sigma_v^2);
    
    % population IV strength
    
    IV_strength(i_spec) = variance_explained_IV / variance_projection_error;
    
end

end

