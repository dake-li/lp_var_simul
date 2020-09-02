function [IRF,nlags,lambda_opt] = LP_shrink_est(data_sim,settings);

% preparations

IRF_hor    = settings.est.IRF_hor;
n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
est_n_lag_BIC  = settings.est.est_n_lag_BIC;
n_lags_fix = settings.est.n_lags_fix;

response_pos = settings.est.IRF_response_var_pos;

with_shock = settings.est.with_shock;
recursive_shock = settings.est.recursive_shock;
with_IV = settings.est.with_IV;

if recursive_shock == 1
    recursive_shock_pos = settings.est.recursive_shock_pos;
end
if with_IV == 1
    IV_est_normalize_var_pos = settings.est.IV_est_normalize_var_pos;
end

lambda = settings.est.lambda;
lambdaRange = settings.est.lambdaRange;
irfLimitOrder = settings.est.irfLimitOrder;

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

% data for LP routine
[y, x, w, H_min, H_max, r] = LP_Penalize_data(Y,recursiveShock,responseV,nlags,IRF_hor - 1,irfLimitOrder);

% leave-one-out cross validation
nT = size(Y,1);
lambdaRange = lambdaRange * nT;
lambdaRange = [1e-4, lambdaRange, 1e10]; % allow regular OLS or completely smoothed
rss_cv = locproj_cv(y, x, w, H_min, H_max, r, lambdaRange, settings.est.CV_folds);
[~,lambda_opt_loc] = min(rss_cv);
lambda_opt = lambdaRange(lambda_opt_loc);

% re-estimate
IRF_resp = locproj(y, x, w, H_min, H_max, r, lambda_opt);
IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0);
IRF = IRF_resp / IRF_normalize;

end