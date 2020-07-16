function [data_sim_new] = residualize_z_fn(data_sim,settings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
est_n_lag_BIC  = settings.est.est_n_lag_BIC;
n_lags_fix = settings.est.n_lags_fix;

% collect data

Y = data_sim.data_y;
Z = data_sim.data_z;
H = [Z, Y];

% set lag length

if est_n_lag == 0
    nlags = n_lags_fix;
elseif est_n_lag_BIC == 1
    [BIC,~] = IC_VAR(H,n_lags_max);
    [~,nlags] = min(BIC);
else
    [~,AIC] = IC_VAR(H,n_lags_max);
    [~,nlags] = min(AIC);
end

% estimate VAR

[~,~,~,~,H_Res] = VAR(H,nlags);

Z_Res = [zeros(nlags,1);H_Res(:,1)];

data_sim_new.data_z = Z_Res;

data_sim_new.data_y = data_sim.data_y;
data_sim_new.data_shock = data_sim.data_shock;

end

