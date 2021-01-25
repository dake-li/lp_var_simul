%%TO DO: REWRITE THE ENTIRE STRUCTURE
clear all;
addpath('Plotting_Functions');

% Plot bias/stddev conditional on observable diagnostics (statistics)


%% Settings

% select lag length specifications
lags_select    = 1;

% select experiments
exper_select = 2;

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

% Apply shared settings
settings_shared;

% Diagnostics
diagns = {'largest_root', 'svar';
          'LM_stat', 'svar';
          'F_stat', 'svar_iv';
          'Granger_stat', 'svar'}; % Names of diagnostics and their sub-fields
diagns_name = {'largroot', 'lm', 'F', 'Granger'}; % Names of diagnostics for file names
diag_qs = [0.25 0.75]; % Lower and upper quantile for conditioning
title_qs = {'$<$Q1', '$>$Q3'}; % Title of lower and upper quantiles in plot


%% FIGURES

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
        
        % Load results
        load_results;
        
        % write down the index of quantiles across MCs
        qs_idx = 2 + find(ismember(res.settings.simul.quantiles,diag_qs)); % index of median number in the quantile list (including mean and std)

        
        %----------------------------------------------------------------
        % Compute Reporting Results
        %----------------------------------------------------------------
        
        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf  = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        the_irf = extract_struct(res.results.irf);
        the_irf = the_irf(:,:,:,methods_select{ne});
        the_irf_err = the_irf - permute(the_true_irf,[1 3 2]);
        the_irf_err_rel = the_irf_err./permute(the_rms_irf, [1 3 2]);
        
        for d=1:length(diagns)

            if ~isfield(res.results, diagns{d,1})
                continue;
            end

            the_diag = res.results.(diagns{d,1}).(diagns{d,2});
            the_diag_q = quantile(median(the_diag,1),diag_qs); % Quantile (across specs) of median (across repetitions)

            f_bias = figure('Position', [100 100 1000 400]);
            f_std = figure('Position', [100 100 1000 400]);
            ax_bias = [];
            ax_std = [];
            
            for iq=1:2 % For lower and upper quantile...
                
                % Select relevant repetitions from sample distribution of the diagnostic
                switch iq
                    case 1
                        the_select = (the_diag<the_diag_q(1)); % Obs. in lower tail
                    case 2
                        the_select = (the_diag>the_diag_q(2)); % Obs. in upper tail
                end

                % Conditional loss
                the_irf_err_rel_sel = the_irf_err_rel;
                the_irf_err_rel_sel(true(size(the_irf_err_rel))==permute(~the_select,[3 1 2 4]))=nan; % Set irrelevant obs. to NaN
                the_bias_rel_cond = squeeze(abs(nanmean(the_irf_err_rel_sel,2))); % Conditional bias
                the_std_rel_cond = squeeze(nanstd(the_irf_err_rel_sel,1,2)); % Conditional std dev
                
                % Conditional bias
                figure(f_bias);
                subplot(1,2,iq);
                plot_loss(horzs(2:end)-1, squeeze(nanmedian(the_bias_rel_cond(2:end,:,:), 2)), [], ...
                    strjoin({exper_plotname, ': Relative Bias', diagns_name{d}, title_qs{iq}}), methods_names_plot, font_size, true);
                same_ylim(ax_bias,gca); % Enforce same ylim
                ax_bias = gca;

                % Conditional std dev
                figure(f_std);
                subplot(1,2,iq);
                plot_loss(horzs(2:end)-1, squeeze(nanmedian(the_std_rel_cond(2:end,:,:), 2)), [], ...
                    strjoin({exper_plotname, ': Relative Std', diagns_name{d}, title_qs{iq}}), methods_names_plot, font_size, true);
                same_ylim(ax_std,gca); % Enforce same ylim
                ax_std = gca;
                

            end
            
            figure(f_bias);
            plot_save(fullfile(output_folder, strcat('diagn_', diagns_name{d}, '_bias_reltruth')), output_suffix);
            figure(f_std);
            plot_save(fullfile(output_folder, strcat('diagn_', diagns_name{d}, '_std_reltruth')), output_suffix);

        end
        
    end
    
end

function same_ylim(ax1, ax2)
    % Enforce same ylim
    if isempty(ax1)
        return;
    end
    ylim1 = get(ax1, 'YLim');
    ylim2 = get(ax2, 'YLim');
    ylim_new = [min(ylim1(1),ylim2(1)) max(ylim1(2),ylim2(2))];
    set(ax1, 'YLim', ylim_new);
    set(ax2, 'YLim', ylim_new);
end