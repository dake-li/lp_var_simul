%% PLOT SIMULATION RESULTS: BIAS/VARIANCE LOSSES
% this version: 09/22/2020

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')

%% SETTINGS

% select lag length specifications
lags_select    = 2;

% select experiments
exper_select = [2,5];

% number of adjacent experiments to combine
n_exper_combine = 2;

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

% Apply shared settings
settings_shared;


%% FIGURES

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
        
        % index in the combined figure panel
        if mod(ne, n_exper_combine) == 0
            i_exper_combine = n_exper_combine;
        else
            i_exper_combine = mod(ne, n_exper_combine);
        end
        
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
            
            if i_exper_combine == 1 % initialize the combined figure panel
                fig_loss(j) = figure('Position', [100 100 500*n_exper_combine 500]);
                fig_rank(j) = figure('Position', [100 100 500*n_exper_combine 500]);
            end

            the_result = sqrt(extract_struct(res.results.(the_objects{j})));
            the_result = the_result(:,:,methods_select{ne});
            the_ranks = permute(tiedrank(permute(the_result, [3 1 2])), [2 3 1]); % Rank procedures from lowest to highest (break ties by averaging)

            % Loss
            figure(fig_loss(j));
            subplot(1,n_exper_combine,i_exper_combine);
            plot_loss(horzs(2:end)-1, squeeze(median(the_result(2:end,:,:)./the_rms_irf, 2)), [], ...
                strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot, font_size, true);
            ax_loss(j, i_exper_combine) = gca;
            
            if i_exper_combine == n_exper_combine % save figure
                same_ylim(ax_loss(j,:)); % enforce same ylim
                add_legend(ax_loss(j,:), methods_names_plot, font_size); % add legend for panel
                plot_save(fullfile(output_folder, strcat('loss_', lower(the_titles{j}), '_reltruth')), output_suffix);
            end
            
            % Average rank
            figure(fig_rank(j));
            subplot(1,n_exper_combine,i_exper_combine);
            plot_loss(horzs(2:end)-1, squeeze(mean(the_ranks(2:end,:,:), 2)), [], ...
                strjoin({exper_plotname, ': Average rank of', the_titles{j}}), methods_names_plot, font_size, true);
            ax_rank(j, i_exper_combine) = gca;
            
            if i_exper_combine == n_exper_combine
                same_ylim(ax_rank(j,:)); % enforce same ylim
                add_legend(ax_rank(j,:), methods_names_plot, font_size); % add legend for panel
                plot_save(fullfile(output_folder, strcat('loss_', lower(the_titles{j}), '_avgrank')), output_suffix);
            end
            
        end
        
    end
    
end