%% COMBINE RESULTS FOR MULTIPLE DGP CHOICE SETS
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf

clear all;
addpath(genpath(fullfile('..', 'Auxiliary_Functions')))
addpath(genpath('Subroutines'))

%% SET UP DESTINATION FOLDER AND FILES

% Experiment setup

spec_id_array = 1; % specification choice set id array
dgp_type = 'G'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)
mode_type = 1; % robustness check mode:
               % 1 (baseline), 2 (small sample), 3 (large sample),
               % 4 (salient series), 5 (more observables), 6 (first diff)

% Cumulative IRF for certain variables

cum_irf_by_trans_code = [0,0,0,0,0,0]; % cumulative IRF or not for each transformation code
% (1) y = x, (2) y = (1-L)x, (3) y = (1-L)^2 x,
% (4) y = ln(x), (5) y = (1-L)ln(x), (6) y = (1-L)^2 ln(x)

if mode_type == 5 % rewrite settings for cumulative IRF
    cum_irf_by_trans_code = [0,1,1,0,1,1]; % cumulative IRF for some transformation codes
end

% Summary statistics across Monte Carlo simulations

summ_option.winsor_percent = 0.025; % winsorize each tail with this percentage to compute winsorized mean and std
summ_option.quantiles = [0.1,0.25,0.5,0.75,0.9]; % summarize MCs in the order of mean, std, winsorized mean, winsorized std, and these quantiles
summ_option.summ_stat_name = [{'mean','std','winsorized_mean','winsorized_std'},...
    strcat('quant_', strsplit(num2str(summ_option.quantiles)))]; % the name of each summary statistic

% Storage folder for results

save_pre = 'Results';

mode_list   = {'baseline', 'small', 'large', 'salient', 'more', 'diff'};
save_mode_dir = mode_list{mode_type}; % set up directory for robustness-check modes

if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));

%% COMBINE ALL THE RESULTS

% Combine results across all spec_ids

res = [];

for spec_id = spec_id_array
    
    % Load one partition
    
    filename = fullfile(save_folder, strcat(strcat('DFM_', dgp_type, '_', estimand_type), '_', num2str(spec_id)));
    res_part = load(filename);
    
    % Compute cumulative IRF for certain variables 
    
    res_part = compute_cum_irf(res_part, cum_irf_by_trans_code);
    
    % Summarize across MCs and combine with previous partitions
    
    res = combine_struct(res, res_part, summ_option, []);
    
end

% Update settings

res.settings.specifications.spec_id_array = spec_id_array;

%% COMPUTE MSE, BIAS, VCE

% Unpack results

DFM_estimate = res.DFM_estimate;
DF_model = res.DF_model;
settings = res.settings;
results = res.results;

% Compute bias2, variance based on summary statistics

[results.MSE, results.BIAS2, results.VCE] = irf_perform_summary(results.irf, DF_model.target_irf, settings);

%% SAVE ALL THE RESULTS

% Save combined results

save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
    'DFM_estimate','DF_model','settings','results','-v7.3');

clear save_folder save_pre save_mode_dir save_suff filename spec_id* summ_option res res_part dgp_type estimand_type lag_type mode_type mode_list cum_irf_by_trans_code