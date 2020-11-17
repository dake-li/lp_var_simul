function [IRF,nlags,F_stat,F_pvalue] = SVAR_IV_est(data_sim,settings);

% residualize IV (reg Z on lagged Z and Y)

data_sim = residualize_z_fn(data_sim,settings);

% preparations

run('Estimation_Setup');

% estimate VAR

[VARout,IVout] = SVAR_IV(Y,IV_est_normalize_var_pos,nlags);
ShockVector = IVout.gamma;
IRF = IRF_SVAR(VARout.By,ShockVector,IRF_hor - 1);
IRF = IRF(response_pos,:)';
F_stat = IVout.Fstat_z; % Wald stat
F_pvalue = chi2cdf(F_stat, 1, 'upper');

end