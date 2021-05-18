%% DFM SIMULATION STUDY: PLOT BIAS/VARIANCE LOSSES
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')
addpath(genpath(fullfile('..', 'Subroutines')))

%% SETTINGS

% select lag length specifications
lags_select    = 2; % options: 1 (AIC), 2 (4 lags), 3 (8 lags)

% select and group experiments
exper_select_group = {[2,5], [3,6], [1,4]}; % combine G and MP for observed shock, recursive, and IV

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 6 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

% select a subset of DGPs
select_DGP = 0; % if select a subset of DGPs?
select_DGP_fn = @(i_dgp, res) res.DF_model.VAR_largest_root(i_dgp) > median(res.DF_model.VAR_largest_root); % binary selection criteria

% Apply shared settings
settings_shared;

%% FIGURES

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
               
        %----------------------------------------------------------------
        % Load Results
        %----------------------------------------------------------------
        
        load_results;
        
        % see if ready to plot for this group of experiments
        if exper_group_end(ne) == 0
            continue;
        end
        
        % keep only the selected subset of DGPs
        if select_DGP == 1
            DGP_selected = arrayfun(@(x) select_DGP_fn(x,res), 1:res.settings.specifications.n_spec)'; % binary DGP selection label
            res = combine_struct(res,[],[],DGP_selected);
        end
        
        %----------------------------------------------------------------
        % Compute Reporting Results
        %----------------------------------------------------------------
        
        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf  = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        % Compute robust statistics
        q1_idx = stat_index(0.25, res.settings); % Index of first quartile
        med_idx = stat_index(0.5, res.settings); % Index of median
        q3_idx = stat_index(0.75, res.settings); % Index of third quartile
        the_fields = fieldnames(res.results.irf);
        for ii=1:length(the_fields)
            res.results.medBIAS2.(the_fields{ii}) = (squeeze(res.results.irf.(the_fields{ii})(:,med_idx,:))-the_true_irf).^2; % Median bias squared
            res.results.IQR2.(the_fields{ii}) = squeeze(res.results.irf.(the_fields{ii})(:,q3_idx,:)-res.results.irf.(the_fields{ii})(:,q1_idx,:)).^2; % IQR squared
        end
        
        %----------------------------------------------------------------
        % Plot Results
        %----------------------------------------------------------------
        
        the_objects = {'BIAS2',    'VCE',   'medBIAS2',     'IQR2'}; % Objects to plot
        the_titles =  {'Bias',     'Std',   'MedBias',      'IQR'};  % Plot titles/file names

        the_methods_index = cellfun(@(x) find(strcmp(res.settings.est.methods_name, x)), methods_fields{ne}); % index of each method

        for j=1:length(the_objects)
            
            the_result = sqrt(extract_struct(res.results.(the_objects{j})));
            the_result = the_result(:,:,the_methods_index);
            the_ranks = permute(tiedrank(permute(the_result, [3 1 2])), [2 3 1]); % Rank procedures from lowest to highest (break ties by averaging)

            % normalized losses
            
            plot_loss(horzs-1, squeeze(median(the_result./the_rms_irf, 2)), [], ...
                strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot, font_size);
            plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_loss_', lower(the_titles{j}), '_reltruth')), output_suffix);
            
            % loss function ranks
            
            plot_loss(horzs-1, squeeze(mean(the_ranks, 2)), [], ...
                strjoin({exper_plotname, ': Average rank of', the_titles{j}}), methods_names_plot, font_size);
            plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_loss_', lower(the_titles{j}), '_avgrank')), output_suffix);
            
        end
        
    end
    
end