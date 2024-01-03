% origin folder for results

rootfolder = fullfile('..', 'Results'); % Root folder with files

% select robustness check mode

mode_list   = {'baseline', 'small', 'large', 'salient', 'more', 'diff'};
mode_folders = mode_list(mode_select);

% select lag length specifications

lags_list      = {'lag_aic', 'lag4', 'lag8', 'lag12'}; % Folders with files
lags_folders   = lags_list(lags_select);

% record experiment group info

exper_select  = cell2mat(exper_select_group); % unzip the groups of experiments
exper_group_num = cellfun('length', exper_select_group); % count number of experiments in each group
exper_group_end = zeros(1, length(exper_select));
exper_group_end(cumsum(exper_group_num)) = 1; % indicate if this experiment is end of group
if size(exper_select_group{1},2) > 1
    indic_grouped = 1; % indicator for whether we are doing grouping
else
    indic_grouped = 0;
end

% select experiments

exper_list        = {'DFM_G_IV', 'DFM_G_ObsShock', 'DFM_G_Recursive', ...
                        'DFM_MP_IV', 'DFM_MP_ObsShock', 'DFM_MP_Recursive'}; % Files in each of the above folders
                    
if indic_grouped == 0
    folder_list       = exper_list;
    exper_names_list  = {'G IV', 'G ObsShock', 'G Recursive', ...
                        'MP IV', 'MP ObsShock', 'MP Recursive'}; % Experiment names for plots
elseif indic_grouped == 1
    folder_list       = {'DFM_IV', 'DFM_ObsShock', 'DFM_Recursive', ...
                        'DFM_IV', 'DFM_ObsShock', 'DFM_Recursive'};
    exper_names_list  = {'IV', 'ObsShock', 'Recursive', ...
                        'IV', 'ObsShock', 'Recursive'}; % Experiment names for plots
end

exper_files   = exper_list(exper_select);
exper_folders = folder_list(exper_select);
exper_names   = exper_names_list(exper_select);

% select estimation methods for each experiment

methods_obsshock_fields  = {'svar','svar_corrbias','bvar','lp','lp_corrbias','lp_penalize','var_avg'};
methods_recursive_fields = methods_obsshock_fields;
methods_iv_fields        = [methods_obsshock_fields {'svar_iv'}];

methods_obsshock_names  = {'VAR','BC VAR','BVAR','LP','BC LP','Pen LP','VAR Avg'};
methods_recursive_names = methods_obsshock_names;
methods_iv_names        = [methods_obsshock_names {'SVAR-IV'}];

methods_select = {methods_iv_select,methods_obsshock_select,methods_recursive_select,...
                    methods_iv_select,methods_obsshock_select,methods_recursive_select};
methods_select = methods_select(exper_select);

methods_fields = {methods_iv_fields(methods_iv_select),methods_obsshock_fields(methods_obsshock_select),methods_recursive_fields(methods_recursive_select),...
                    methods_iv_fields(methods_iv_select),methods_obsshock_fields(methods_obsshock_select),methods_recursive_fields(methods_recursive_select)};
methods_fields = methods_fields(exper_select);

methods_names = {methods_iv_names(methods_iv_select),methods_obsshock_names(methods_obsshock_select),methods_recursive_names(methods_recursive_select),...
                    methods_iv_names(methods_iv_select),methods_obsshock_names(methods_obsshock_select),methods_recursive_names(methods_recursive_select)};
methods_names = methods_names(exper_select);

% select DGP subsets

if DGP_select == 1
    select_DGP_fn = @(i_dgp, res) any(res.settings.specifications.var_select(i_dgp,:)>=res.settings.specifications.random_category_range(11,1)); % binary selection criteria: specifications with asset price & sentiment
elseif DGP_select == 2
    select_DGP_fn = @(i_dgp, res) res.DF_model.R0_sq(i_dgp) <= prctile(res.DF_model.R0_sq,10); % binary selection criteria: low degree of invertibility
elseif DGP_select == 3
    select_DGP_fn = @(i_dgp, res) res.DF_model.R0_sq(i_dgp) >= prctile(res.DF_model.R0_sq,90); % binary selection criteria: high degree of invertibility
end

%----------------------------------------------------------------
% Figure Output
%----------------------------------------------------------------

font_size = 16;         % Font size of method labels (titles and tick labels will be scaled accordingly)
output_dir = 'fig';     % Folder
output_suffix = 'png';  % File suffix