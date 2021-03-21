function [IRF,nlags,weightOut] = VAR_avg_est(data_sim,settings)
% Function for estimating IRFs using VAR model averaging

% preparations

average_max_lags = settings.est.average_max_lags;
h_store = settings.est.average_store_weight;
run('Estimation_Setup'); % common setup for all estimation methods
weightOut = NaN(2 * n_lags_max, length(h_store)); % placeholder for weights on different VAR submodels

% check if averaging from one lag to max number of lags

if average_max_lags == 1
    n_lags_average_over = n_lags_max; % to max number of lags
else
    n_lags_average_over = nlags; % to fixed/estimated number of lags
end

% estimate different VAR submodels

[combination_irf,weights,G] = VAR_ModelAverage(Y,recursiveShock,responseV,n_lags_average_over,IRF_hor - 1,settings.est.average_options);
IRF = combination_irf / G(normalizeV, recursiveShock); % IRF normalize by response of normalization variable. Warning: correspond to matrix C in our paper

% store weights

weightOut(1:(n_lags_average_over * 2), :) = weights(:, h_store - 1); % record weights of submodels at horizon h

end

