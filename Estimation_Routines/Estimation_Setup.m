%% Common Setup for All Estimation Methods

% unpack settings

IRF_hor    = settings.est.IRF_hor;
n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
est_n_lag_BIC  = settings.est.est_n_lag_BIC;
n_lags_fix = settings.est.n_lags_fix;

response_pos = settings.est.IRF_response_var_pos;
normalize_pos = settings.est.est_normalize_var_pos;

with_shock = settings.est.with_shock;
recursive_shock = settings.est.recursive_shock;
with_IV = settings.est.with_IV;

if recursive_shock == 1
    recursive_shock_pos = settings.est.recursive_shock_pos;
end

% collect data

if with_shock == 1 % observe shock: w_t = (shock, \bar{w}_t)
    Y = [data_sim.data_shock,data_sim.data_y]; % Warning: correspond to w_t in our paper
    responseV = response_pos + 1; % location of response variable
    recursiveShock = 1; % location of impulse variable
    normalizeV = normalize_pos + 1; % location of normalization variable
elseif with_IV == 1 % IV: w_t = (IV, \bar{w}_t)
    Y = [data_sim.data_z,data_sim.data_y];
    responseV = response_pos + 1;
    recursiveShock = 1;
    normalizeV = normalize_pos + 1;
else % recursive: w_t = \bar{w}_t
    Y = data_sim.data_y;
    responseV = response_pos;
    recursiveShock = recursive_shock_pos;
    normalizeV = normalize_pos;
end

% set lag length

if est_n_lag == 0 % fix lag order
    nlags = n_lags_fix;
elseif est_n_lag_BIC == 1 % estimate lag order via BIC
    [BIC,~] = IC_VAR(Y,n_lags_max);
    [~,nlags] = min(BIC);
else % estimate lag order via AIC
    [~,AIC] = IC_VAR(Y,n_lags_max);
    [~,nlags] = min(AIC);
end