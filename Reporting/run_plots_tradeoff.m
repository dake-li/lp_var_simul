%% PLOT SIMULATION RESULTS: BIAS/VARIANCE TRADE-OFF
% this version: 09/22/2020

%% HOUSEKEEPING

clc
clear all
close all

addpath('Auxiliary_Functions')
warning('off','MATLAB:structOnObject')

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

%----------------------------------------------------------------
% Preference Plots
%----------------------------------------------------------------

% bias weight grid

n_weight    = 10001;
weight_grid = linspace(1,0,n_weight)';

% reference method (should have something here on this being LP...)

base_name = {'LP'};
base_indic = NaN(length(exper_files),1);

for ne=1:length(exper_files)
    if sum(strcmp(methods_names{ne}, base_name)) == 0
        error('The base method is not included!')
    else
        base_indic(ne) = find(strcmp(methods_names{ne}, base_name));
    end
end

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

n = 200;
clear cmap
cmap(1,:) = [1 1 1];
cmap(2,:) = [0.5 0.5 0.5];
cmap(3,:) = [0 0 0];

[X,Y] = meshgrid([1:3],[1:50]);

cmap = interp2(X([1,25,50],:),Y([1,25,50],:),cmap,X,Y);

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
        base_method_name = methods_names_plot{base_indic(ne)};
        
        %----------------------------------------------------------------
        % Compute Reporting Results
        %----------------------------------------------------------------
        
        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        the_objects = methods_names_plot; % Objects to plot
        the_titles =  methods_names_plot; % Plot titles/file names

        for j=1:length(the_objects)
            
            if j ~= base_indic(ne)

            the_BIAS2 = extract_struct(res.results.BIAS2);
            the_BIAS2 = the_BIAS2(:,:,methods_select{ne});
            the_VCE   = extract_struct(res.results.VCE);
            the_VCE   = the_VCE(:,:,methods_select{ne});
            
            the_BIAS2rel = the_BIAS2./the_rms_irf.^2;
            the_VCErel   = the_VCE./the_rms_irf.^2;
            
            pref_base = zeros(n_weight,max(res.settings.est.IRF_select),size(the_BIAS2rel,2));
            for i_weight = 1:n_weight
                loss_base   = weight_grid(i_weight) * the_BIAS2rel(:,:,base_indic(ne)) + (1-weight_grid(i_weight)) * the_VCErel(:,:,base_indic(ne));
                loss_method = weight_grid(i_weight) * the_BIAS2rel(:,:,j) + (1-weight_grid(i_weight)) * the_VCErel(:,:,j);

                pref_base(i_weight,:,:) = (loss_base <= loss_method);
            end
            pref_base = mean(pref_base,3);

            plot_tradeoff(pref_base(:,2:end), cmap, horzs(2:end)-1, weight_grid, ...
                strjoin({exper_plotname, ': Trade-Off for', the_titles{j}, 'vs', base_method_name}))
            plot_save(fullfile(output_folder, strcat(lower(the_titles{j}), '_tradeoff')), output_suffix);

            end

        end
        
    end
    
end