%% COMBINE RESULTS FOR MULTIPLE SPECIFICATION CHOICE SETS


%% SET UP DESTINATION FOLDER AND FILES

save_pre = fullfile('..', 'Results');

spec_id_array = 1:2; % specification choice set id array
dgp_type = 'G'; % 'MP'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % 'Recursive'; 'IV'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)

if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, strcat('lag', save_suff));


%% COMBINE ALL THE RESULTS

% combine results across all spec_ids

[DFM_estimate, DF_model, settings, results] = combine_struct(save_folder, strcat('DFM_', dgp_type, '_', estimand_type), spec_id_array);

% save combined results

save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
    'DFM_estimate','DF_model','settings','results',...
    'spec_id_array','dgp_type','estimand_type','lag_type','-v7.3');

clear save_folder save_pre save_suff