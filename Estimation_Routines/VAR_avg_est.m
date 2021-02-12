function [IRF,nlags,weightOut] = VAR_avg_est(data_sim,settings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% preparations

average_max_lags = settings.est.average_max_lags;
h_store = settings.est.average_store_weight;
run('Estimation_Setup');
weightOut = NaN(2 * n_lags_max, length(h_store));

% check if averaging over max number of models

if average_max_lags == 1
    n_lags_average_over = n_lags_max;
else
    n_lags_average_over = nlags;
end

% estimate VAR

[combination_irf,weights,G] = VAR_ModelAverage(Y,recursiveShock,responseV,n_lags_average_over,IRF_hor - 1,settings.est.average_options);
IRF = combination_irf / G(normalizeV, recursiveShock);

% store weights
weightOut(1:(n_lags_average_over * 2), :) = weights(:, h_store - 1); % record weights of submodels at h

end

