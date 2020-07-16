function [LRV_Cov_tr_ratio, VAR_largest_root, frac_coef_for_large_lags] = compute_persist_DSGE(model, settings)
%   compute CovMat: 
%       using ABCDEF representation
%   compute Innovation Representation:
%       following state-space form (A,B,C,D) based on Fernandez et al. (2005)
%   compute LRVMat:
%       use VMA in innovation representation and tranform spectral density
%   compute LRV_Cov_tr_ratio:
%       ratio of trace(LRV) over trace(Cov) for each specification
%   compute truncated VAR representation:
%       truncate infinite order VAR in innovation representation at a large lag order
%   compute VAR_largest_root:
%       largest root using truncated VAR for each specification
%   compute frac_coef_for_large_lags:
%       ratio btw VAR coefficient summed up from lag p+1 to lag infinity over 
%       VAR coefficient summed up from lag 1 to lag infinity using
%       truncated VAR for each specification

% prepare

var_select = settings.specifications.var_select;
n_spec = size(var_select,1);
n_var = size(var_select,2);

VAR_infinity_truncate = settings.est.VAR_infinity_truncate; % truncate infinite-order VAR
VAR_fit_nlags = settings.est.n_lags_fix; % examine population fit for VAR(p)

LRV_Cov_tr_ratio = NaN(n_spec, 1);
VAR_largest_root = NaN(n_spec, 1);
frac_coef_for_large_lags = NaN(n_spec,1);

for i_spec = 1:n_spec
    
    A   = model.ABCD.A;
    B   = model.ABCD.B;
    C   = model.ABCD.C;
    C   = C(var_select(i_spec,:),:);
    D   = model.ABCD.D;
    D   = D(var_select(i_spec,:),:);
    n_s = size(A, 1);
    n_y = size(C, 1);
   
    %----------------------------------------------------------------
    % Compute Covariance Matrix for All Observables
    %----------------------------------------------------------------

    % covariance matrix for s

    BB = B * B';
    vec_BB = BB(:);
    vec_CovMat_s = (eye(n_s^2) - kron(A, A)) \ vec_BB;
    CovMat_s = reshape(vec_CovMat_s, [n_s, n_s]);

    % covariance matrix for y
    CovMat_y = C * CovMat_s * C' + D * D';

    %----------------------------------------------------------------
    % Compute Innovation Representation
    %----------------------------------------------------------------

    % derive innovation representation
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

    % innovation for y

    Sigma_u = C * Sigma_states * C' + D * D';
    
    % compute cholesky decomposition
    
    G = chol(Sigma_u, 'lower');
    
    %----------------------------------------------------------------
    % Compute LRV Matrix for All Observables
    %----------------------------------------------------------------

    % LRV of white noise

    LRVMat_WN = eye(n_y);

    % transformation function Theta
    
    Theta = (G + C * ((eye(n_s) - A) \ K_states) * G);

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
    
    % compute VAR polynomial for y
    VAR_poly = NaN(n_var, n_var, 1+VAR_infinity_truncate);
    VAR_poly(:,:,1) = eye(n_var);
    for ilag = 1:VAR_infinity_truncate
        VAR_poly(:,:,1+ilag) = - C * (A - K_states*C)^(ilag-1) * K_states;
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