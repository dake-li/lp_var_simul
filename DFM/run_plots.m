clear all;

% Plot simulation results


%% Settings

% Results .mat files
mat_rootfolder = 'Results'; % Root folder with files
mat_folders = {'lag4', 'lag8'}; % Folders with files
mat_files = {'DFM_G_IV', 'DFM_G_ObsShock', 'DFM_G_Recursive', ...
             'DFM_MP_IV', 'DFM_MP_ObsShock', 'DFM_MP_Recursive'}; % Files in each of the above folders

% Output settings for figures
output_dir = 'Figures'; % Folder
output_suffix = 'png';  % File suffix

% Plotting
plot_smooth_param = 0.1; % Smoothing parameter for cubic spline


%% Create plots for each experiment

for nf=1:length(mat_folders) % For each folder...

    for ne=1:length(mat_files) % For each experiment in folder...
        
        %% File/folder names
        
        exper_name = mat_files{ne}; % Name of current experiment
        file_name = fullfile(mat_folders{nf}, exper_name); % Name of .mat results file
        output_folder = fullfile(output_dir, file_name); % Name of output folder
        mkdir(output_folder); % Create output folder
        

        %% Load results

        res = load(fullfile(mat_rootfolder, strcat(file_name, '.mat'))); % Load
        horzs = res.settings.est.IRF_select; % Impulse response horizons
        methods_name = res.settings.est.methods_name; % Names of estimation methods


        %% Report features of DGP

        % Degree of invertibility
        figure;
        histogram(res.DF_model.R0_sq, 'Normalization', 'probability');
        title(strjoin({exper_name, ': degree of invertibility'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'R0sq'), output_suffix);
        
        % Persistence
        figure;
%         histogram(mean(res.DF_model.persistency,2), 'Normalization', 'probability');
%         title(strjoin({exper_name, ': average persistence'}), 'Interpreter', 'none');
        histogram(res.DF_model.LRV_Cov_tr_ratio, 'Normalization', 'probability');
        title(strjoin({exper_name, ': LRV to Var ratio'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'persistence'), output_suffix);
        
        % IV strength
        if isfield(res.DF_model, 'IV_strength')
            figure;
            histogram(res.DF_model.IV_strength, 'Normalization', 'probability');
            title(strjoin({exper_name, ': IV strength'}), 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'IVstrength'), output_suffix);
        end
        
        
        %% Estimated tuning parameters

        % Number of lags
        the_nlags = res.results.n_lags.svar;
        figure;
        histogram(the_nlags(:), 'Normalization', 'probability');
        title(strjoin({exper_name, ': number of lags (across specs+sims)'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'nlags'), output_suffix);
        
        figure;
        histogram(std(the_nlags), 'Normalization', 'probability');
        title(strjoin({exper_name, ': std (across sims) of number of lags'}), 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'nlags_std'), output_suffix);

%         % Shrinkage penalty
%         the_lambda = res.results.lambda.lp_penalize;
%         figure;
%         histogram(the_lambda(:), 'Normalization', 'probability');
%         title(strjoin({exper_name, ': shrinkage penalty (across specs+sims)'}), 'Interpreter', 'none');
%         plot_save(fullfile(output_folder, 'lambda'), output_suffix);
%         
%         figure;
%         histogram(std(the_lambda), 'Normalization', 'probability');
%         title(strjoin({exper_name, ': std (across sims) of shrinkage penalty'}), 'Interpreter', 'none');
%         plot_save(fullfile(output_folder, 'lambda_std'), output_suffix);
        
        % Model-averaging weights
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
        sgtitle(strjoin({exper_name, ': average model-avg weight (across specs+sims)'}), 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
        plot_save(fullfile(output_folder, 'weight'), output_suffix);


        %% MSE/bias/var

        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        the_objects = {'MSE',   'BIAS2',    'VCE'}; % Objects to plot
        the_titles =  {'RMSE',  'bias',     'std'}; % Plot titles/file names

        for j=1:length(the_objects)

            the_result = sqrt(extract_struct(res.results.(the_objects{j}))); % Extract results for object of interest
            the_ranks = permute(tiedrank(permute(the_result, [3 1 2])), [2 3 1]); % Rank procedures from lowest to highest (break ties by averaging)

            % Relative to true IRF
            plot_results(horzs, squeeze(median(the_result./the_rms_irf, 2)), [], strjoin({exper_name, ': median of', the_titles{j}, 'relative to true IRF'}), methods_name);
            plot_save(fullfile(output_folder, strcat(lower(the_titles{j}), '_reltruth')), output_suffix);

            % Relative to VAR
            plot_results(horzs, squeeze(median(the_result(:,:,2:end)./the_result(:,:,1), 2)), 1, strjoin({exper_name, ': median of', the_titles{j}, 'relative to VAR'}), methods_name(2:end));
            plot_save(fullfile(output_folder, strcat(lower(the_titles{j}), '_relvar')), output_suffix);

            % Average rank
            plot_results(horzs, squeeze(mean(the_ranks, 2)), [], strjoin({exper_name, ': average rank of', the_titles{j}}), methods_name);
            plot_save(fullfile(output_folder, strcat(lower(the_titles{j}), '_avgrank')), output_suffix);

            % Fraction smallest
            plot_results(horzs, squeeze(mean(the_ranks==1, 2)), [], strjoin({exper_name, ': fraction smallest', the_titles{j}}), methods_name);
            plot_save(fullfile(output_folder, strcat(lower(the_titles{j}), '_fracsmallest')), output_suffix);

        end
        
        
        %% MSE vs. persistence
        
%         the_mse = extract_struct(res.results.MSE);
%         plot_smooth(mean(res.DF_model.persistency,2), squeeze(sqrt(mean(the_mse))./the_rms_irf), plot_smooth_param, strjoin({exper_name, ': RMSE vs. persistence (across specs)'}), methods_name, 'average persistence', 'RMSE across horizons relative to true IRF');
%         plot_save(fullfile(output_folder, 'rmse_vspersistence'), output_suffix);
        
        
        %% IV F-stat
        
        if isfield(res.results, 'F_stat')
            
            the_Fstats = res.results.F_stat.svar_iv;
            the_rse_acrosshorz = sqrt(squeeze(mean((extract_struct(res.results.irf)-permute(the_true_irf, [1 3 2])).^2,1))); % Root squared estimation error averaged across horizons (for each spec+sim)
            
            % Average F stat
            figure;
            histogram(the_Fstats(:), 'Normalization', 'probability');
            title(strjoin({exper_name, ': F stat (across specs+sims)'}), 'Interpreter', 'none');
            plot_save(fullfile(output_folder, 'F'), output_suffix);
            
            % RSE vs. F stat scatter
            plot_smooth(the_Fstats(:), reshape(the_rse_acrosshorz./the_rms_irf, [], length(methods_name)), plot_smooth_param, strjoin({exper_name, ': RSE vs. F stat (across specs+sims)'}), methods_name, 'F stat', 'RSE across horizons relative to true IRF');
            plot_save(fullfile(output_folder, 'F_vsrse'), output_suffix);
            
        end


    end

end


%% Auxiliary functions

function [M, fields] = extract_struct(the_struct)

    % Function for extracting numerical results from struct

    fields = fieldnames(the_struct);
    dim = size(the_struct.(fields{1}));
    M = nan([dim length(fields)]);
    for i=1:length(fields)
        if length(dim) == 2
            M(:,:,i) = the_struct.(fields{i});
        elseif length(dim) == 3
            M(:,:,:,i) = the_struct.(fields{i});
        end
    end

end

function plot_results(horzs, results, add_line, plot_name, plot_legend)

    % Function for plotting results across estimation methods

    line_colors = repmat(lines(5),2,1);
    line_styles = [repmat({'-'},5,1); repmat({'--'},5,1)];
    
    figure;
    
    hold on;
    for i=1:size(results,2)
        plot(horzs, results(:,i)', 'Color', line_colors(i,:), 'LineStyle', line_styles{i});
    end
    hold off;
    
    if ~isempty(add_line)
        hold on;
        the_xlim = xlim;
        line(the_xlim, add_line*ones(1,2), 'Color', 'k', 'LineStyle', ':');
        hold off;
        xlim(the_xlim);
    end
    title(plot_name, 'Interpreter', 'none');
    legend(plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'Interpreter', 'none');

end

function plot_smooth(x, y, smooth_param, plot_name, plot_legend, plot_xlab, plot_ylab)

    xeval = linspace(quantile(x,0.01),quantile(x,0.99),5e3)';
    smooth_y = zeros(5e3,size(y,2));
    for j=1:size(y,2)
        smooth_y(:,j) = csaps(x, y(:,j), smooth_param, xeval);
    end
    
    figure;
    plot(xeval, smooth_y);
    title(plot_name, 'Interpreter', 'none');
    xlabel(plot_xlab);
    ylabel(plot_ylab);
    legend(plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'Interpreter', 'none');

end

function plot_save(filename, suffix)

    % Function for saving figures
    
    saveas(gcf, strcat(filename, '.', suffix));
    close(gcf);

end
