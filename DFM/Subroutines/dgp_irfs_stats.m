function model = dgp_irfs_stats(model, settings, estimand_type)

% Impulse responses and summary statistics for DGPs


var_select = settings.specifications.var_select;
[n_spec,n_var] = size(var_select);

model.LRV_Cov_tr_ratio = nan(n_spec,1);
model.dLRV_dCov_tr_ratio = nan(n_spec,1);
model.VAR_largest_root = nan(n_spec,1);
model.VAR_quant_root = nan(n_spec,1);
model.frac_coef_for_large_lags = nan(n_spec,1);
model.R0_sq = nan(n_spec,1);
model.IV_strength = nan(n_spec,1);

%% True IRFs in encompassing DFM

model.irf = compute_irfs(model.ABCD,settings.est.shock_weight,settings.est.IRF_hor);

if any(strcmp(estimand_type, {'ObsShock', 'IV'})) % Observed shock and IV
    model.normalized_irf = compute_normalized_irfs(model,settings);
    model.target_irf = model.normalized_irf(settings.est.IRF_select, :);
else % Recursive
    model.target_irf = nan(settings.est.IRF_hor,n_spec); % To be computed below
end

% Loop over DGPs
for i_spec = 1:n_spec

    %% ABCD model for given subset of observables
    ABCD_obs = model.ABCD;
    ABCD_obs.C = ABCD_obs.C(var_select(i_spec,:),:);
    ABCD_obs.D = ABCD_obs.D(var_select(i_spec,:),:);

    %% Compute reduced-form VAR(infinity) representation

    % Reduced-dimensionality ABCD model for selected observables
    [ABCD_small, shock_select] = ABCD_reduce(ABCD_obs);

    % VAR(infinity)
    red_form = reduced_form_VAR(ABCD_small,settings.est.VAR_infinity_truncate);

    % Cholesky impulse responses for recursive identification
    if strcmp(estimand_type, 'Recursive')
        G = chol(red_form.innov_var, 'lower');
        chol_irfs = compute_irfs(red_form.innov_ABCD,G(:,settings.est.recursive_shock_pos),settings.est.IRF_hor);
        model.target_irf(:,i_spec) = chol_irfs(:,settings.est.IRF_response_var_pos)/chol_irfs(1,settings.est.est_normalize_var_pos);
    end

    %% Reduced-form summary statistics

    % Roots of reduced-form VAR
    nlag_comp = settings.est.VAR_infinity_truncate_comp; % # lags to include in companion matrix
    comp_form = [cell2mat(red_form.coef(1:nlag_comp)); eye(n_var*(nlag_comp-1)) zeros(n_var*(nlag_comp-1), n_var)]; % Companion matrix
    roots = abs(eig(comp_form));
    model.VAR_largest_root(i_spec) = max(roots); % Max
    model.VAR_quant_root(i_spec) = quantile(roots,settings.est.VAR_root_quant); % Percentile

    % tr(LRV)/tr(Var) ratio
    if model.VAR_largest_root(i_spec)<0.999
        [cov, lrv] = cov_lrv(ABCD_small);
        model.LRV_Cov_tr_ratio(i_spec) = trace(lrv)/trace(cov);
    else
        dABCD = ABCD_diff(ABCD_small); % ABCD representation of first differences
        [dcov, dlrv] = cov_lrv(dABCD);
        model.dLRV_dCov_tr_ratio(i_spec) = trace(dlrv)/trace(dcov);
    end

    % Fraction of VAR coefficients at long lags
    norms = cellfun(@(x) norm(x,'fro'), red_form.coef); % Frobenius norm of each VAR coefficient matrix at lags 1,2,...
    model.frac_coef_for_large_lags(i_spec) = 1-sum(norms(1:settings.est.n_lags_fix))/sum(norms);

    %% Shock and IV summary statistics

    % Degree of invertibility
    if ~strcmp(estimand_type, 'Recursive')
        model.R0_sq(i_spec) = degree_invertibility(ABCD_small.D, red_form.innov_var, settings.est.shock_weight(shock_select));
    else
        model.R0_sq(i_spec) = 1;
    end

    % IV strength
    if strcmp(estimand_type, 'IV')
        model.IV_strength(i_spec) = IV_strength(ABCD_obs, model.IV, settings);
    end

end


end