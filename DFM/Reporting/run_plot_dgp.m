%% DFM SIMULATION STUDY: TABLES AND PLOTS FOR DGP FEATURES
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clc
close all
clear all;
addpath('Plotting_Functions');
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

% Summary statistics for table
tab_stat = {'LRV_Cov_tr_ratio', 'VAR_largest_root', 'frac_coef_for_large_lags', 'R0_sq'}; % Summary stats to copy (in addition to IRF stats defined below)
tab_quants = [0.1 0.25 0.5 0.75 0.9]; % Quantiles to report across specifications

% Table with model specification tests
spec_lags_cutoff = 2; % Report how often quantiles of lag length exceeds this number
spec_lags_quants = [0.5 0.75 0.9]; % Quantiles of lag length to report
spec_lm_power = 0.25; % Report fraction of DGPs with at least this power of LM test
spec_lm_signifs = [0.05 0.1 0.25]; % Significance levels for LM test

% IRF examples to plot
spec_select = [10 20 30 3010 3020 3030];
linestyles = {'-', '--', ':', '-o', '--o', ':o'};
colors = lines(6);

%% CREATE TABLES AND PLOTS

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
        
        % Record the index of median across MCs
        median_idx = stat_index(0.5, res.settings); % index of median number in the quantile list (including mean and std)        
        
        %----------------------------------------------------------------
        % Tables of Summary Statistics
        %----------------------------------------------------------------
        
        tab = table;
        
        % Basic DGP summary statistics
        for is=1:length(tab_stat)
            tab.(tab_stat{is}) = res.DF_model.(tab_stat{is});
        end
        
        % IV F-stat
        if isfield(res.DF_model, 'IV_strength')
            tab.iv_F = res.settings.simul.T./(1./res.DF_model.IV_strength-1);
        end
        
        % IRF summary statistics
        tab.irf_num_local_extrema = sum(diff(sign(diff(res.DF_model.target_irf',1,2)),1,2)~=0,2);
        [~,I] = max(abs(res.DF_model.target_irf)',[],2);
        tab.irf_maxabs_h = res.settings.est.IRF_select(I)'-1;
        tab.irf_mean_max_ratio = mean(res.DF_model.target_irf',2)./max(abs(res.DF_model.target_irf)',[],2);
        
        % R-squared from regressing IRFs on quadratic
        X = [ones(res.settings.est.IRF_hor,1) (res.settings.est.IRF_select').^(1:2)];
        betas = X\res.DF_model.target_irf;
        resid = res.DF_model.target_irf-X*betas;
        tab.irf_R2s = 1-(var(resid)./var(res.DF_model.target_irf))';
        
        % Create table of quantiles
        tab_summ = varfun(@(x) [min(x) quantile(x,tab_quants) max(x)]', tab);
        tab_summ.Properties.VariableNames=regexprep(tab_summ.Properties.VariableNames, 'Fun_', '');
        tab_summ2 = table;
        tab_summ2.quantile = [0; tab_quants(:); 1];
        tab_summ = [tab_summ2 tab_summ];
        
        % Write to file
        writetable(tab_summ, fullfile(output_folder, 'dgp_summ.csv'));
        
        clearvars I X betas resid tab tab_summ tab_summ2;        
        
        %----------------------------------------------------------------
        % Tables of Specification Tests
        %----------------------------------------------------------------
        
        tab_spec = table;

        % Number of lags
        for ii=1:length(spec_lags_quants)
            tab_spec.(sprintf('%s%02d', 'nlag_exceed_q', round(100*spec_lags_quants(ii)))) ...
                = mean(res.results.n_lags.svar(stat_index(spec_lags_quants(ii), res.settings),:)>=spec_lags_cutoff);
        end
        
        % Power of LM test
        for ii=1:length(spec_lm_signifs)
            tab_spec.(sprintf('%s%02d', 'lm_power_', round(100*spec_lm_signifs(ii)))) ...
                = mean(res.results.LM_pvalue.svar(stat_index(spec_lm_power, res.settings),:)<spec_lm_signifs(ii));
        end
        
        % Write to file
        writetable(tab_spec, fullfile(output_folder, 'dgp_spec.csv'));
        
        clearvars tab_spec;
        
        
        %----------------------------------------------------------------
        % Plots to summarize DGP Features
        %----------------------------------------------------------------
                
        % Degree of invertibility
        figure;
        histogram(res.DF_model.R0_sq, 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': degree of invertibility'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_R0sq'), output_suffix);
        
        % Persistence
        figure;
        histogram(res.DF_model.LRV_Cov_tr_ratio, 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': LRV to Var ratio'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_LRVVar'), output_suffix);
        
        figure;
        histogram(res.DF_model.VAR_largest_root, 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': largest VAR root'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_largroot'), output_suffix);
        
        figure;
        histogram(res.DF_model.frac_coef_for_large_lags, 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': fraction of long-lag VAR coefs'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_longlag'), output_suffix);
        
        % IV strength
        if isfield(res.DF_model, 'IV_strength')
            figure;
            histogram(res.DF_model.IV_strength, 'Normalization', 'probability');
            title(strjoin({exper_plotname, ': IV strength'}), 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'dgp_IVstrength'), output_suffix);
        end
        
        % Some IRFs
        figure('Units', 'inches', 'Position', [0 0 8 4]);
        hold on;
        norm_irf = @(x) x/max(abs(x));
        for i_spec_indx = 1:length(spec_select)
            i_spec = spec_select(i_spec_indx);
            plot(res.settings.est.IRF_select-1, norm_irf(res.DF_model.target_irf(:,i_spec)), ...
                 linestyles{i_spec_indx}, 'Color', colors(i_spec_indx,:), 'Linewidth', 2, ...
                 'MarkerSize', 4, 'MarkerFaceColor', colors(i_spec_indx,:));
        end
        hold off;
        xlabel('Horizon','interpreter','latex','FontSize',12);
        set(gca,'XTick',[min(res.settings.est.IRF_select-1) 2:2:max(res.settings.est.IRF_select-1)]);
        xlim([min(res.settings.est.IRF_select-1) max(res.settings.est.IRF_select-1)]);
        grid on;
        set(gca,'TickLabelInterpreter','latex');
        set(gca,'FontSize',12);
        plot_save(fullfile(output_folder, 'dgp_irfs'), output_suffix);        
        
        %----------------------------------------------------------------
        % Tuning Parameters for Estimation Methods
        %----------------------------------------------------------------
        
        % Number of lags
        the_nlags = res.results.n_lags.svar;
        figure;
        histogram(the_nlags(median_idx,:), 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': median number of lags (across specs)'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_nlags'), output_suffix);
        
        figure;
        histogram(the_nlags(2,:), 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': std (across sims) of number of lags'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_nlags_std'), output_suffix);

        % Shrinkage penalty
        the_lambda = res.results.lambda.lp_penalize;
        figure;
        [~,the_edges] = histcounts(log10(the_lambda(median_idx,:)));
        histogram(the_lambda(median_idx,:),10.^the_edges);
        set(gca, 'xscale','log'); % Log scale for x axis
        title(strjoin({exper_plotname, ': median shrinkage penalty (across specs)'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_lambda'), output_suffix);
        
        % Model-averaging weights
        if isfield(res.results, 'weight')
            
            the_weights = res.results.weight.var_avg;
            the_store_weights = res.settings.est.average_store_weight;
            the_maxlag = res.settings.est.n_lags_max;
            
            figure;
            
            for j=1:length(the_store_weights) % For each horizon where weights are stored...
                subplot(1,length(the_store_weights),j);
                plot(1:the_maxlag, mean(reshape(the_weights(1:the_maxlag,j,1,:), the_maxlag, []), 2)); % AR weights
                hold on;
                plot(1:the_maxlag, mean(reshape(the_weights(the_maxlag+1:end,j,1,:), the_maxlag, []), 2)); % VAR weights
                hold off;
                title(sprintf('%s%d', 'h = ', the_store_weights(j)));
                xlabel('no. of lags');
                legend('AR', 'VAR', 'Location', 'NorthEast');
            end
            sgtitle(strjoin({exper_plotname, ': average model-avg weight (across specs+sims)'}), 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'dgp_weight'), output_suffix);
            
        end       
        
        %----------------------------------------------------------------
        % IV F-Statistics
        %----------------------------------------------------------------
        
        if isfield(res.results, 'F_stat')
                       
            the_Fstats = res.results.F_stat.svar_iv;
            
            % Average F stat
            figure;
            histogram(the_Fstats(1,:), 'Normalization', 'probability');
            title(strjoin({exper_plotname, ': average F stat (across specs)'}), 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'dgp_F'), output_suffix);
            
        end
        
        %----------------------------------------------------------------
        % Granger Causality
        %----------------------------------------------------------------
        
        if isfield(res.results, 'Granger_stat')
            
            the_Grangerstats = res.results.Granger_stat.svar;
            
            % Average Granger causality Wald-stat
            figure;
            histogram(the_Grangerstats(1,:), 'Normalization', 'probability');
            title(strjoin({exper_plotname, ': average Granger causality Wald stat (across specs)'}), 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'dgp_Granger'), output_suffix);
            
        end


    end

end