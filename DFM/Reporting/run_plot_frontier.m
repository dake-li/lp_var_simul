%% DFM SIMULATION STUDY: PLOT BIAS/VARIANCE FRONTIER
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')
addpath(genpath(fullfile('..', 'Subroutines')))

%% SETTINGS

% select robustness check mode
mode_select    = 1; % options: 1 (baseline), 2 (cumulative IRF), 3 (persistent DGP), 4 (small sample), 5 (salient series)

% select lag length specifications
lags_select    = 2; % options: 1 (AIC), 2 (4 lags), 3 (8 lags)

% select and group experiments
exper_select_group = {[2,5]}; % combine G and MP for observed shock, recursive, and IV

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 6 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

% select a subset of DGPs
select_DGP = 0; % if select a subset of DGPs?
select_DGP_fn = @(i_dgp, res) res.DF_model.VAR_largest_root(i_dgp) > median(res.DF_model.VAR_largest_root); % binary selection criteria

% report quantile loss across DGPs
loss_quant = 0.5; % report which quantile loss across DGPs? (default is median loss, i.e. 0.5)

% Apply shared settings
settings_shared;

% Plot settings
horzs_plot = [4 10]; % Plot results for these horizons side by side
method_labels = {'VAR', [0.06 0.04], [0.06 0.04];
                 'Bias-corrected VAR', [0.1 -0.2], [-0.04 -0.09];
                 'Bayesian VAR', [0 0.07], [0 0.07];
                 'LP', [0 0.05], [0 0.05];
                 'Penalized LP', [0.05 0], [0.05 0];
                 'VAR model averaging', [0 0.05], [0 0.05]};
    % Label for each estimator, as well as coordinate offset (in
    % 'normalized' units) for text box relative to data point
fontsize_plot = 12; % Font size for plot axes


%% FIGURES

for n_mode=1:length(mode_folders) % For each robustness check mode...

    for nf=1:length(lags_folders) % For each lag-order folder...
    
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
            
            % bias and standard deviation
            
            the_bias = squeeze(quantile(sqrt(extract_struct(res.results.BIAS2))./the_rms_irf, loss_quant, 2));
            the_std = squeeze(quantile(sqrt(extract_struct(res.results.VCE))./the_rms_irf, loss_quant, 2));  
            
            horzs_sel = ismember(res.settings.est.IRF_select, horzs_plot); % Selected horizons
            the_bias_sel = the_bias(horzs_sel,:);
            the_std_sel = the_std(horzs_sel,:);
            
            % approximate frontier for visual aid
            
            the_bias_frontier_grid = linspace(0.8 * min(the_bias_sel(:)),1.2 * max(the_bias_sel(:)),100)';
            the_max_std = max(the_std_sel(:))*1.2;
            
            %----------------------------------------------------------------
            % Plot Results
            %----------------------------------------------------------------
            
            if loss_quant == 0.5
                remark_loss_quant = ''; % remark in file name for quantile loss
            else
                remark_loss_quant = strcat('_p', num2str(round(loss_quant*100)));
            end

            for ih=1:length(horzs_plot)
                figure('Units', 'normalize', 'Position', [ih*0.3 0.3 0.3 0.5]);
                plot_frontier(the_bias_sel(ih,:)',the_std_sel(ih,:)',the_bias_frontier_grid,the_max_std,...
                    [method_labels(:,1) method_labels(:,ih+1)],fontsize_plot,...
                    sprintf('%s%d%s', 'Horizon: ', horzs_plot(ih), ' Quarters'));
                set(gcf, 'PaperPositionMode', 'auto');
                plot_save(fullfile(output_folder, sprintf('%s%d', 'frontier_h', horzs_plot(ih), remark_loss_quant)), output_suffix);
            end
            
        end
        
    end

end