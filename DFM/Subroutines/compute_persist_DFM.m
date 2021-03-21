function [LRV_Cov_tr_ratio, VAR_largest_root, frac_coef_for_large_lags] = compute_persist_DFM(model, settings)
% Function for computing different population measures of persistence of \bar{w}_t in each DGP
    %-------------------------------------------------
    % Long-run variance over variance trace-ratio:
    %-------------------------------------------------
    % first, directly use ABCDEF representation to compute Var(\bar{w}_t):
        % ABCDEF representation:
            % state transition:  s_t = A * s_{t-1} + B * \epsilon_t
            % measurement eq:    y_t = C * s_{t-1} + D * \epsilon_t + e^*_t
            % measurement error: e_t = E * e_{t-1} + F * \omega_t

            % Warning: y_t corresponds to X_t in the paper

    % second, transform DFM representation into ABCD representation:  
        % ABCD representation: (See details in Technical Companion Note)
            % state transition:  s_{t+1} = A * s_t + B * \zeta_t
            % measurement eq:    \bar{w}^*_t = C * s_t + D * \zeta_t

            % where \bar{w}^*_t = (1 - \Delta(L)) * \bar{w}_t

    % third, transform ABCD representation into Wold representation:
        % first, VAR(\infty) representation:
            % \bar{w}^*_t = \Psi^*(L) \bar{w}^*_{t-1} + u_t (See details in Fernandez et al., 2005)
        % second, VAR(\infty) representation:
            % \bar{w}_t = \Psi(L) \bar{w}_{t-1} + u_t
        % then, VMA(infty) representation:
            % \bar{w}_t = \Theta(L) u_t

    % fourth, compute LRV(\bar{w}_t) = \Theta(1) * LRV(u_t) * \Theta(1)

    % finally, compute tr(LRV) / tr(Var) for each DGP
        
    %-------------------------------------------------
    % VAR largest root:
    %-------------------------------------------------
    % re-use VAR(\infty) representation:
        % \bar{w}_t = \Psi(L) \bar{w}_{t-1} + u_t
        
    % truncate infinite order VAR in innovation representation at a large lag order
            
    % compute largest root using truncated VAR for each DGP
                
    %-------------------------------------------------
    % Fraction of VAR coef. beyond lag p:
    %-------------------------------------------------
    % re-use VAR(\infty) representation:
        % \bar{w}_t = \Psi(L) \bar{w}_{t-1} + u_t
        
    % sum up VAR coefficients from lag p+1 to lag infinity
    % (approximated by truncated VAR at a large lag order)
    
    % sum up VAR coefficients from lag 1 to lag infinity
    % (approximated by truncated VAR at a large lag order)
    
    % compute the ratio of these two terms for each DGP

%-------------------------------------------------
% Prepare
%-------------------------------------------------

% unpack settings

n_y = model.n_y;
n_s = model.n_s;
model_A = model.ABCD.A;
model_B = model.ABCD.B;
model_C = model.ABCD.C;
model_D = model.ABCD.D;
model_E = model.ABCD.E;
model_F = model.ABCD.F;

n_fac = model.n_fac;
n_lags_fac = model.n_lags_fac;
n_lags_uar = model.n_lags_uar;
n_lags_state = max(n_lags_uar + 1, n_lags_fac);

Phi = model.Phi;
Sigma_eta = model.Sigma_eta;
Lambda = model.Lambda;
delta = model.delta;
sigma_v = model.sigma_v;

var_select = settings.specifications.var_select;
n_spec = size(var_select,1);
n_var = size(var_select,2);

VAR_infinity_truncate = settings.est.VAR_infinity_truncate; % truncate infinite-order VAR at this lag length
VAR_fit_nlags = settings.est.n_lags_fix; % later examine population fit for VAR(p)

LRV_Cov_tr_ratio = NaN(n_spec, 1);
VAR_largest_root = NaN(n_spec, 1);
frac_coef_for_large_lags = NaN(n_spec,1);

% go thru each DGP

for i_spec = 1:n_spec
   
    %----------------------------------------------------------------
    % Compute Variance-Covariance Matrix for All Observables
    %----------------------------------------------------------------

    % covariance matrix for s

    model_BB = model_B * model_B';
    vec_model_BB = model_BB(:);
    vec_CovMat_s = (eye(n_s^2) - kron(model_A, model_A)) \ vec_model_BB;
    CovMat_s = reshape(vec_CovMat_s, [n_s, n_s]);

    % covariance matrix for e
    
    EF_row_index = zeros(n_var, n_lags_uar);
    for i_lag_uar = 1:n_lags_uar
        EF_row_index(:, i_lag_uar) = var_select(i_spec,:) + (i_lag_uar - 1) * n_y;
    end
    EF_row_index = EF_row_index(:);
    model_FF = model_F(EF_row_index, var_select(i_spec,:)) * model_F(EF_row_index, var_select(i_spec,:))';
    vec_model_FF = model_FF(:);
    vec_CovMat_e = (eye((n_var * n_lags_uar)^2) - ...
        kron(model_E(EF_row_index, EF_row_index), model_E(EF_row_index, EF_row_index))) \ vec_model_FF;
    CovMat_e = reshape(vec_CovMat_e, [n_var * n_lags_uar, n_var * n_lags_uar]);
    CovMat_e_star = CovMat_e(1:n_var, 1:n_var);

    % covariance matrix for y
    CovMat_y = model_C(var_select(i_spec,:), :) * CovMat_s * model_C(var_select(i_spec,:), :)' +...
        model_D(var_select(i_spec,:), :) * model_D(var_select(i_spec,:), :)' + CovMat_e_star;

    %----------------------------------------------------------------
    % Compute ABCD Representation
    %----------------------------------------------------------------

    % compute A

    A = zeros(n_lags_state * n_fac);
    A(1:n_fac, 1:(n_lags_fac * n_fac)) = Phi(1:n_fac, :);
    A((1 + n_fac):(n_lags_state * n_fac), 1:((n_lags_state - 1) * n_fac)) = eye((n_lags_state - 1) * n_fac);

    % compute B

    B = zeros(n_lags_state * n_fac, n_fac + n_var);
    B(1:n_fac, 1:n_fac) = chol(Sigma_eta, 'lower');

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
    
    %----------------------------------------------------------------
    % Compute VAR(\infty) Representation
    %----------------------------------------------------------------

    % compute steady state conditional variance in Kalman filter

    cond_var = cond_var_fn_St_1(A, B, C, D);

    % compute cov-var matrix of innovations

    Sigma_innovation = C * cond_var * C' + D * D';

    % compute Kalman gain

    K = (A * cond_var * C' + B * D') / Sigma_innovation;

    % compute cholesky decomposition

    G = chol(Sigma_innovation, 'lower');
    
    %----------------------------------------------------------------
    % Compute LRV Matrix for All Observables
    %----------------------------------------------------------------

    % LRV of white noise

    LRVMat_WN = eye(n_var);

    % transformation function Theta

    Theta_left = eye(n_var);
    for ilag = 1:n_lags_uar
        Theta_left = Theta_left - diag(delta(var_select(i_spec,:), ilag));
    end
    Theta_right = (G + C * ((eye(n_lags_state * n_fac) - A) \ K) * G);
    Theta = Theta_left \ Theta_right;

    % LRV of observables

    LRVMat_y = Theta * LRVMat_WN * Theta';

    %----------------------------------------------------------------
    % Compute Ratio between tr(LRV) over tr(Cov)
    %----------------------------------------------------------------

    % compute ratio of tr(LRV) over tr(Cov) in each specification

    LRV_Cov_tr_ratio(i_spec) = trace(LRVMat_y) / trace(CovMat_y);
    
    %----------------------------------------------------------------
    % Compute Truncated VAR in Innovation Representation
    %----------------------------------------------------------------
    
    % compute VAR polynomial for y_star
    VAR_poly_star = NaN(n_var, n_var, 1+VAR_infinity_truncate);
    VAR_poly_star(:,:,1) = eye(n_var);
    for ilag = 1:VAR_infinity_truncate
        VAR_poly_star(:,:,1+ilag) = - C * (A - K*C)^(ilag-1) * K;
    end
    
    % compute VAR polynomial for y
    VAR_poly = VAR_poly_star;
    for ilag_uar = 1:n_lags_uar
        for ilag = 1:(VAR_infinity_truncate+1-ilag_uar)
            VAR_poly(:,:,ilag+ilag_uar) = VAR_poly(:,:,ilag+ilag_uar) +...
                VAR_poly_star(:,:,ilag) * (-diag(delta(var_select(i_spec, :), ilag_uar)));
        end
    end
    
    %----------------------------------------------------------------
    % Compute Largest Root
    %----------------------------------------------------------------
    
    % compute companion-form VAR for y
    VAR_companion_form = zeros(n_var * VAR_infinity_truncate);
    VAR_companion_form(1:n_var,:) = reshape(-VAR_poly(:,:,2:end), [n_var, n_var*VAR_infinity_truncate]);
    VAR_companion_form((n_var+1):end, 1:(n_var*(VAR_infinity_truncate-1))) = eye(n_var*(VAR_infinity_truncate-1));
    
    % compute largest root in VAR
    VAR_largest_root(i_spec) = max(abs(eig(VAR_companion_form)));
    
    %----------------------------------------------------------------
    % Compute VAR(p) Fit by Examining Coefficients after Lag p
    %----------------------------------------------------------------
    
    sum_coef_all_lags = 0;
    sum_coef_for_large_lags = 0;
    for ilag = 1:VAR_infinity_truncate
        sum_coef_all_lags = sum_coef_all_lags + norm(VAR_poly(:,:,ilag),'fro');
        if ilag > VAR_fit_nlags
            sum_coef_for_large_lags = sum_coef_for_large_lags + norm(VAR_poly(:,:,ilag),'fro');
        end
    end
    frac_coef_for_large_lags(i_spec) = sum_coef_for_large_lags / sum_coef_all_lags;
    
end

end