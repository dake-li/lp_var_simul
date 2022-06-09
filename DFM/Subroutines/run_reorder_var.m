%% REORDER VARIABLE LIST IN OLD RESULT FILES TO MATCH S&W 2016
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clear all;
clc;

%% SET UP DESTINATION FOLDER AND FILES

% Experiment setup

dgp_type = 'G'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)

% Reordering setup

reorder    = [1:76, 87:94, 77:86, 95:171, 181:195, 172:180, 196:207]; % index to reorder data to match variable list in Stock-Watson (2016)

% Storage folder for results

input_pre = fullfile('..', 'Results');
output_pre = fullfile('..', 'Results_reorder');
save_mode_dir = 'baseline';
if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
input_folder = fullfile(input_pre, save_mode_dir, strcat('lag', save_suff));
output_folder = fullfile(output_pre, save_mode_dir, strcat('lag', save_suff));

%% REORDER THE RESULTS

% Load result
    
filename = fullfile(input_folder, strcat('DFM_', dgp_type, '_', estimand_type));
res = load(filename);

% Unpack

DFM_estimate = res.DFM_estimate;
DF_model = res.DF_model;
settings = res.settings;
results = res.results;

% Reorder DFM_estimate

DFM_estimate.Lambda(reorder, :) = DFM_estimate.Lambda;
DFM_estimate.sigma_v(reorder) = DFM_estimate.sigma_v;
DFM_estimate.delta(reorder, :) = DFM_estimate.delta;

DFM_estimate.bpnamevec(reorder) = DFM_estimate.bpnamevec;
DFM_estimate.bplabvec_long(reorder) = DFM_estimate.bplabvec_long;
DFM_estimate.bplabvec_short(reorder) = DFM_estimate.bplabvec_short;
DFM_estimate.bptcodevec(reorder) = DFM_estimate.bptcodevec;

DFM_estimate.r2(reorder) = DFM_estimate.r2;

% Reorder DF_model

DF_model.Lambda        = DFM_estimate.Lambda;
DF_model.delta         = DFM_estimate.delta;
DF_model.sigma_v       = DFM_estimate.sigma_v;

DF_model.variable_name = DFM_estimate.bplabvec_long;
DF_model.trans_code = DFM_estimate.bptcodevec;

DF_model.ABCD  = ABCD_fun_DFM(DF_model);

DF_model.irf(:, reorder) = DF_model.irf;

% Reorder settings

settings.specifications.var_select = reshape(reorder(settings.specifications.var_select(:)), ...
    [settings.specifications.n_spec, settings.specifications.n_var]);

% Store reordering info

DF_model.reorder = reorder;

%% SAVE ALL THE RESULTS

% Save reordered results

mkdir(output_folder);
save(fullfile(output_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
    'DFM_estimate','DF_model','settings','results','-v7.3');