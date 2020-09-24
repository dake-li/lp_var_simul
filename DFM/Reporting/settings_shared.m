% file origin

rootfolder = fullfile('..', 'Results'); % Root folder with files

% select lag length specifications

lags_list      = {'lag4','lag8'}; % Folders with files
lags_folders   = lags_list(lags_select);

% select experiments

exper_list        = {'DFM_G_IV', 'DFM_G_ObsShock', 'DFM_G_Recursive', ...
                        'DFM_MP_IV', 'DFM_MP_ObsShock', 'DFM_MP_Recursive'}; % Files in each of the above folders
exper_names_list  = {'G IV', 'G ObsShock', 'G Recursive', ...
                        'MP IV', 'MP ObsShock', 'MP Recursive'}; % Experiment names for plots

exper_files   = exper_list(exper_select);
exper_names   = exper_names_list(exper_select);

% select estimation methods for each experiment

methods_iv_names        = {'SVAR','Bias-Corr. SVAR','BVAR','LP','Penalized LP','VAR Avg.','SVAR-IV'};
methods_obsshock_names  = {'SVAR','Bias-Corr. SVAR','BVAR','LP','Penalized LP','VAR Avg.'};
methods_recursive_names = {'SVAR','Bias-Corr. SVAR','BVAR','LP','Penalized LP','VAR Avg.'};

methods_select = {methods_iv_select,methods_obsshock_select,methods_recursive_select,...
                    methods_iv_select,methods_obsshock_select,methods_recursive_select};
methods_select = methods_select(exper_select);

methods_names = {methods_iv_names(methods_iv_select),methods_obsshock_names(methods_obsshock_select),methods_recursive_names(methods_recursive_select),...
                    methods_iv_names(methods_iv_select),methods_obsshock_names(methods_obsshock_select),methods_recursive_names(methods_recursive_select)};
methods_names = methods_names(exper_select);

%----------------------------------------------------------------
% Figure Output
%----------------------------------------------------------------

output_dir = 'fig';     % Folder
output_suffix = 'png';  % File suffix