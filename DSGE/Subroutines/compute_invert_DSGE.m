function R0_sq = compute_invert_DSGE(model,settings);

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

var_select = settings.specifications.var_select;
n_spec = settings.specifications.n_spec;

shock_weight = settings.est.shock_weight;

R0_sq = NaN(n_spec, 1);

for i_spec = 1:n_spec
    
    A   = model.ABCD.A;
    B   = model.ABCD.B;
    C   = model.ABCD.C;
    C   = C(var_select(i_spec,:),:);
    D   = model.ABCD.D;
    D   = D(var_select(i_spec,:),:);
    n_s = size(A, 1);
    n_w = size(B, 2);
    n_y = size(C, 1);

    %----------------------------------------------------------------
    % Link Reduced-Form Errors and Structural Shocks
    %----------------------------------------------------------------

    % derive Sigma and K

    Sigma_states = eye(n_s);
    K_states = zeros(n_s,n_y);

    dist = 1;
    tol = 10^(-10);
    relax = 0.9;

    while dist >= tol
        Sigma_upd = (A - K_states*C)*Sigma_states*(A-K_states*C)' + B*B' + K_states*(D*D')*K_states' - B * D'*K_states' - K_states * D * B';
        K_upd = (A * Sigma_upd * C' + B * D') * (C * Sigma_upd * C' + D * D')^(-1);

        dist1 = max(max(abs(Sigma_upd - Sigma_states)));
        dist2 = max(max(abs(K_upd - K_states)));

        Sigma_states = relax * Sigma_states + (1-relax) * Sigma_upd;
        K_states = relax * K_states + (1-relax) * K_upd;

        dist = max(dist1, dist2);
    end

    % get the M matrices

    M = NaN(n_y, n_w, 1);
    M(:,:,1) = D;

    %----------------------------------------------------------------
    % Total Available Weights
    %----------------------------------------------------------------

    dg_invert_cov_mat = M(:,:,1)' * (C * Sigma_states * C' + D * D')^(-1) * M(:,:,1);

    %----------------------------------------------------------------
    % Final R_0^2 for true structural shock
    %----------------------------------------------------------------

    R0_sq(i_spec) = shock_weight' * dg_invert_cov_mat * shock_weight;

end

end