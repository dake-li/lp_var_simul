%% PLOT SIMULATION RESULTS: BIAS/VARIANCE LOSSES
% this version: 09/22/2020

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')

%% SETTINGS

% select lag length specifications
lags_select    = [1 2];

% select experiments
exper_select = [2 6];

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

% Apply shared settings
settings_shared;


%% FIGURES

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
        
        % Load results
        load_results;
        
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
            the_ranks = permute(tiedrank(permute(the_result, [3 1 2])), [2 3 1]); % Rank procedures from lowest to highest (break ties by averaging)

            % Loss
            plot_loss(horzs(2:end)-1, squeeze(median(the_result(2:end,:,:)./the_rms_irf, 2)), [], ...
                strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot, font_size);
            plot_save(fullfile(output_folder, strcat('loss_', lower(the_titles{j}), '_reltruth')), output_suffix);
            
            % Average rank
            plot_loss(horzs(2:end)-1, squeeze(mean(the_ranks(2:end,:,:), 2)), [], ...
                strjoin({exper_plotname, ': Average rank of', the_titles{j}}), methods_names_plot, font_size);
            plot_save(fullfile(output_folder, strcat('loss_', lower(the_titles{j}), '_avgrank')), output_suffix);

        end
        
    end
    
end