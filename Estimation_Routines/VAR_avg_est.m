function [IRF,nlags,weightOut] = VAR_avg_est(data_sim,settings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% preparations

IRF_hor    = settings.est.IRF_hor;
n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
est_n_lag_BIC  = settings.est.est_n_lag_BIC;
n_lags_fix = settings.est.n_lags_fix;

average_max_lags = settings.est.average_max_lags;
h_store = settings.est.average_store_weight;

response_pos = settings.est.IRF_response_var_pos;

with_shock = settings.est.with_shock;
recursive_shock = settings.est.recursive_shock;
with_IV = settings.est.with_IV;

weightOut = NaN(2 * n_lags_max, length(h_store));

if recursive_shock == 1
    recursive_shock_pos = settings.est.recursive_shock_pos;
end
if with_IV == 1
    IV_est_normalize_var_pos = settings.est.IV_est_normalize_var_pos;
end

% collect data

if with_shock == 1
    Y = [data_sim.data_shock,data_sim.data_y];
    responseV = response_pos + 1;
    recursiveShock = 1;
    normalizeV = recursiveShock;
elseif with_IV == 1
    Y = [data_sim.data_z,data_sim.data_y];
    responseV = response_pos + 1;
    recursiveShock = 1;
    normalizeV = IV_est_normalize_var_pos + 1;
else
    Y = data_sim.data_y;
    responseV = response_pos;
    recursiveShock = recursive_shock_pos;
    normalizeV = recursiveShock;
end

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

% check if averaging over max number of models

if average_max_lags == 1
    n_lags_average_over = n_lags_max;
else
    n_lags_average_over = nlags;
end

% estimate VAR

[combination_irf,weights,G] = VAR_ModelAverage(Y,recursiveShock,responseV,n_lags_average_over,IRF_hor - 1);
IRF = combination_irf / G(normalizeV, recursiveShock);

% store weights
weightOut(1:(n_lags_average_over * 2), :) = weights(:, h_store - 1); % record weights of submodels at h

end

