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

% select robustness check mode
mode_select    = 1; % options: 1 (baseline), 2 (cumulative IRF), 3 (persistent DGP), 4 (small sample), 5 (salient series)

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

% regress bias/std on variable category counts
reg_cat = 1; % if run regression?
reg_cat_horz = []; % if non-empty, only use subset of horizons for regression (e.g., [1 2] means first and second estimated horizons)

% report quantile loss across DGPs
loss_quant = 0.5; % report which quantile loss across DGPs? (default is median loss, i.e. 0.5)

% Apply shared settings
settings_shared;

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
            % Prepare regression on variable category counts (if desired)
            %----------------------------------------------------------------
            
            if reg_cat==1
                
                % Mapping between variables and categories
                aux = zeros(res.DF_model.n_y,1);
                aux(res.settings.specifications.random_category_range(:,1)) = 1;
                var_cat = cumsum(aux); % mapping between variables and categories
                spec_cat = var_cat(res.settings.specifications.var_select); % list of categories included in each DGP
                
                % Categories in each DGP
                num_cat = size(res.settings.specifications.random_category_range,1); % no. of categories
                aux2 = reshape(bsxfun(@eq,spec_cat(:),1:num_cat),res.settings.specifications.n_spec,res.settings.specifications.n_var,num_cat);
                spec_cat_num = permute(sum(aux2,2),[1 3 2]); % category counts for each DGP
                
                % Covariate matrix for regressions
                aux3 = repmat(eye(res.settings.est.IRF_hor),res.settings.specifications.n_spec,1); % indicators for horizon
                the_reg_cat_horz = reg_cat_horz;
                if isempty(the_reg_cat_horz)
                    the_reg_cat_horz = 1:res.settings.est.IRF_hor; % all horizons
                end
                reg_sel = any(aux3(:,the_reg_cat_horz),2); % include only selected horizons in regressions
                aux4 = kron(spec_cat_num,ones(res.settings.est.IRF_hor,1));
                reg_cat_X = [aux4(:,1:end-1) aux3(:,the_reg_cat_horz)]; % omit last category and any undesired horizons
                reg_cat_vars = [strcat('cat', cellfun(@num2str, num2cell(1:num_cat-1), 'UniformOutput', false)) strcat('h', cellfun(@num2str, num2cell(res.settings.est.IRF_select(the_reg_cat_horz)-1), 'UniformOutput', false))];
                reg_cat_vars = reg_cat_vars(:);
                
                clearvars aux aux2 aux3 aux4;
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

            if loss_quant == 0.5
                remark_loss_quant = ''; % remark in file name for quantile loss
            else
                remark_loss_quant = strcat('_p', num2str(round(loss_quant*100)));
            end
    
            the_methods_index = cellfun(@(x) find(strcmp(res.settings.est.methods_name, x)), methods_fields{ne}); % index of each method
    
            for j=1:length(the_objects)
                
                the_result = sqrt(extract_struct(res.results.(the_objects{j})));
                the_result = the_result(:,:,the_methods_index);
                the_ranks = permute(tiedrank(permute(the_result, [3 1 2])), [2 3 1]); % Rank procedures from lowest to highest (break ties by averaging)
    
                % normalized losses
                
                plot_loss(horzs-1, squeeze(quantile(the_result./the_rms_irf, loss_quant, 2)), [], ...
                    strjoin({exper_plotname, ': Relative', the_titles{j}}), methods_names_plot, font_size);
                plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_loss_', lower(the_titles{j}), '_reltruth', remark_loss_quant)), output_suffix);
                
                % loss function ranks
                
                plot_loss(horzs-1, squeeze(mean(the_ranks, 2)), [], ...
                    strjoin({exper_plotname, ': Average rank of', the_titles{j}}), methods_names_plot, font_size);
                plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_loss_', lower(the_titles{j}), '_avgrank')), output_suffix);
                
                % regression on variable category counts
                
                if reg_cat==1
                    reg_cat_Y = log(reshape(the_result./the_rms_irf,res.settings.est.IRF_hor*res.settings.specifications.n_spec,[])); % log loss
                    reg_beta = reg_cat_X(reg_sel,:)\reg_cat_Y(reg_sel,:); % OLS regression of log loss on category and horizon variables
                    the_tab = array2table(reg_beta);
                    the_tab.Properties.VariableNames = res.settings.est.methods_name(the_methods_index);
                    the_tab_var = table;
                    the_tab_var.VARIABLE = reg_cat_vars;
                    writetable([the_tab_var the_tab],fullfile(output_folder, strcat(exper_names{ne}, '_loss_', lower(the_titles{j}), '_regcat', '.csv'))); % write table to file
                end
                
            end
            
        end
        
    end

end