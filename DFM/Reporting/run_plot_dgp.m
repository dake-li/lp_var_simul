clear all;
addpath('Plotting_Functions');

% Plot features of GDP and estimated tuning parameters


%% Settings

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


%% Create plots for each experiment

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
        
        % Load results
        load_results;


        %% Report features of DGP

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
        
        
        %% Estimated tuning parameters

        % Number of lags
        the_nlags = res.results.n_lags.svar;
        figure;
        histogram(the_nlags(:), 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': number of lags (across specs+sims)'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_nlags'), output_suffix);
        
        figure;
        histogram(std(the_nlags), 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': std (across sims) of number of lags'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_nlags_std'), output_suffix);

        % Shrinkage penalty
        the_lambda = res.results.lambda.lp_penalize;
        figure;
        [~,the_edges] = histcounts(log10(the_lambda));
        histogram(the_lambda,10.^the_edges);
        set(gca, 'xscale','log'); % Log scale for x axis
        title(strjoin({exper_plotname, ': shrinkage penalty (across specs+sims)'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_lambda'), output_suffix);
        
        figure;
        histogram(std(log(the_lambda)), 'Normalization', 'probability');
        title(strjoin({exper_plotname, ': std (across sims) of log shrinkage penalty'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'dgp_lambda_logstd'), output_suffix);
        
        % Model-averaging weights
        if isfield(res.results, 'weight')
            
            the_weights = res.results.weight.var_avg;
            the_store_weights = res.settings.est.average_store_weight;
            the_maxlag = res.settings.est.n_lags_max;
            figure;
            for j=1:length(the_store_weights) % For each horizon where weights are stored...
                subplot(1,length(the_store_weights),j);
                plot(1:the_maxlag, mean(reshape(the_weights(1:the_maxlag,j,:,:), the_maxlag, []), 2)); % AR weights
                hold on;
                plot(1:the_maxlag, mean(reshape(the_weights(the_maxlag+1:end,j,:,:), the_maxlag, []), 2)); % VAR weights
                hold off;
                title(sprintf('%s%d', 'h = ', the_store_weights(j)));
                xlabel('no. of lags');
                legend('AR', 'VAR', 'Location', 'NorthEast');
            end
            sgtitle(strjoin({exper_plotname, ': average model-avg weight (across specs+sims)'}), 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'dgp_weight'), output_suffix);
            
        end
        
        
        %% IV F-stat
        
        if isfield(res.results, 'F_stat')
            
            the_Fstats = res.results.F_stat.svar_iv;
            
            % Average F stat
            figure;
            histogram(the_Fstats(:), 'Normalization', 'probability');
            title(strjoin({exper_plotname, ': F stat (across specs+sims)'}), 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'dgp_F'), output_suffix);
            
        end


    end

end
