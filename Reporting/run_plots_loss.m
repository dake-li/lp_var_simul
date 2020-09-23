%% PLOT SIMULATION RESULTS: BIAS/VARIANCE LOSSES
% this version: 09/22/2020

%% HOUSEKEEPING

clc
clear all
close all

addpath('Auxiliary_Functions')

%% SETTINGS

%----------------------------------------------------------------
% Root Folders
%----------------------------------------------------------------

% file origin

rootfolder = '/Users/christianwolf/Dropbox/Research/lp_var_simul/stored_result/small_run_0916'; % Root folder with files

% select lag length specifications

lags_list      = {'lag4','lag8'}; % Folders with files
lags_select    = [1 2];
lags_folders   = lags_list(lags_select);

% select experiments

exper_list        = {'DFM_G_IV', 'DFM_G_ObsShock', 'DFM_G_Recursive', ...
                        'DFM_MP_IV', 'DFM_MP_ObsShock', 'DFM_MP_Recursive'}; % Files in each of the above folders
exper_names_list  = {'G IV', 'G ObsShock', 'G Recursive', ...
                        'MP IV', 'MP ObsShock', 'MP Recursive'}; % Experiment names for plots
         
exper_select = [1 2 3 4 5 6];

exper_files   = exper_list(exper_select);
exper_names   = exper_names_list(exper_select);

% select estimation methods for each experiment

methods_iv_names        = {'SVAR','Bias-Corr. SVAR','BVAR','LP','Penalized LP','VAR Avg.','SVAR-IV'};
methods_obsshock_names  = {'SVAR','Bias-Corr. SVAR','BVAR','LP','Penalized LP','VAR Avg.'};
methods_recursive_names = {'SVAR','Bias-Corr. SVAR','BVAR','LP','Penalized LP','VAR Avg.'};

methods_iv_select        = [1 2 3 4 5 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

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

%% FIGURES

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
        
        %----------------------------------------------------------------
        % File/Folder Names
        %----------------------------------------------------------------
        
        exper_filename = exper_files{ne}; % Name of current experiment
        exper_plotname = exper_names{ne};
        file_name = fullfile(lags_folders{nf}, exper_filename); % Name of .mat results file
        output_folder = fullfile(output_dir, file_name); % Name of output folder
        mkdir(output_folder); % Create output folder        

        %----------------------------------------------------------------
        % Load Simulation Results
        %----------------------------------------------------------------

        res = load(fullfile(rootfolder, strcat(file_name, '.mat'))); % Load
        horzs = res.settings.est.IRF_select; % Impulse response horizons
        methods_names_plot = methods_names{ne};
        
        %----------------------------------------------------------------
        % Compute Reporting Results
        %----------------------------------------------------------------
        
        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf  = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        the_objects = {'MSE',   'BIAS2',    'VCE'}; % Objects to plot
        the_titles =  {'RMSE',  'Bias',     'Std'}; % Plot titles/file names

        for j=1:length(the_objects)

            the_result = sqrt(extract_struct(res.results.(the_objects{j})));
            the_result = the_result(:,:,methods_select{ne});

            plot_loss(horzs(2:end)-1, squeeze(median(the_result(2:end,:,:)./the_rms_irf, 2)), [], ...
                strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot);
            plot_save(fullfile(output_folder, strcat(lower(the_titles{j}), '_reltruth')), output_suffix);

        end
        
    end
    
end