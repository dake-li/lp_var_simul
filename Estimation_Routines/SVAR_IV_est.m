function [IRF,nlags,Fstat] = SVAR_IV_est(data_sim,settings);

% preparations

IRF_hor    = settings.est.IRF_hor;
n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
est_n_lag_BIC  = settings.est.est_n_lag_BIC;
n_lags_fix = settings.est.n_lags_fix;

response_pos = settings.est.IRF_response_var_pos;
normalize_pos = settings.est.IV_est_normalize_var_pos;

% residualize IV (reg Z on lagged Z and Y)

data_sim = residualize_z_fn(data_sim,settings);

% collect data

Y = data_sim.data_y;
Z = data_sim.data_z;
H = [Z,Y]; % combine Z and Y

% set lag length

if est_n_lag == 0
    nlags = n_lags_fix;
elseif est_n_lag_BIC == 1
    [BIC,~] = IC_VAR(Y,n_lags_max);
    [~,nlags] = min(BIC);
else
    [~,AIC] = IC_VAR(Y,n_lags_max);
    [~,nlags] = min(AIC);
end

% estimate VAR

[VARout,IVout] = SVAR_IV(H,normalize_pos,nlags);
ShockVector = IVout.gamma;
IRF = IRF_SVAR(VARout.By,ShockVector,IRF_hor - 1);
IRF = IRF(response_pos,:)';
Fstat = IVout.Fstat_z;

end