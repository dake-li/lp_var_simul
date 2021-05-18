%% COMBINE RESULTS FOR MULTIPLE DGP CHOICE SETS
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clear all;

%% SET UP DESTINATION FOLDER AND FILES

spec_id_array = 1; % specification choice set id array
dgp_type = 'G'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)
summ_option.winsor_percent = 0.025; % winsorize each tail with this percentage to compute winsorized mean and std
summ_option.quantiles = [0.1,0.25,0.5,0.75,0.9]; % summarize MCs in the order of mean, std, winsorized mean, winsorized std, and these quantiles
summ_option.summ_stat_name = [{'mean','std','winsorized_mean','winsorized_std'},strcat('quant_', strsplit(num2str(summ_option.quantiles)))]; % the name of each summary statistic

save_pre = fullfile('..', 'Results');
if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, strcat('lag', save_suff));


%% COMBINE ALL THE RESULTS

% combine results across all spec_ids

res = [];

for spec_id = spec_id_array
    
    filename = fullfile(save_folder, strcat(strcat('DFM_', dgp_type, '_', estimand_type), '_', num2str(spec_id)));
    res_part = load(filename);
    res = combine_struct(res, res_part, summ_option, []);
    
end

res.settings.specifications.spec_id_array = spec_id_array;

%% SAVE ALL THE RESULTS

% unpack results

DFM_estimate = res.DFM_estimate;
DF_model = res.DF_model;
settings = res.settings;
results = res.results;

% save combined results

save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
    'DFM_estimate','DF_model','settings','results','-v7.3');

clear save_folder save_pre save_suff filename spec_id* summ_option res res_part dgp_type estimand_type lag_type