%% COMBINE RESULTS FOR MULTIPLE SPECIFICATION CHOICE SETS
clear all;

%% SET UP DESTINATION FOLDER AND FILES

spec_id_array = 1:2; % specification choice set id array
dgp_type = 'G'; % Either 'G' or 'MP'
estimand_type = 'IV'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = NaN; % No. of lags to impose in estimation, or NaN (meaning AIC)
quantiles = [0.1,0.25,0.5,0.75,0.9]; % summarize MCs in the order of mean, var, and these quantiles

save_pre = fullfile('..', 'Results');
if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, strcat('lag', save_suff));


%% COMBINE ALL THE RESULTS

% combine results across all spec_ids

[DFM_estimate, DF_model, settings, results] = combine_struct(save_folder, strcat('DFM_', dgp_type, '_', estimand_type), spec_id_array, quantiles);

% save combined results

save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
    'DFM_estimate','DF_model','settings','results',...
    'dgp_type','estimand_type','lag_type','-v7.3');

clear save_folder save_pre save_suff