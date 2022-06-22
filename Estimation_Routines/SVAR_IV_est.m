function [IRF,n_lags_est,F_stat,F_pvalue] = SVAR_IV_est(data_sim,settings);
% Function for estimating IRFs using SVAR-IV

% residualize IV (reg Z on lagged Z and Y)

data_sim = residualize_z_fn(data_sim,settings);

% preparations

run('Estimation_Setup'); % common setup for all estimation methods

% estimate SVAR via IV

[VARout,IVout] = SVAR_IV(Y,normalize_pos,nlags);
ShockVector = IVout.gamma;

% estimate IRF

IRF = IRF_SVAR(VARout.By,ShockVector,IRF_hor - 1);
IRF = IRF(response_pos,:)'; % IRF to one unit of shock (already normalized)
F_stat = IVout.Fstat_z; % Wald stat to test IV strength
F_pvalue = chi2cdf(F_stat, 1, 'upper');

end