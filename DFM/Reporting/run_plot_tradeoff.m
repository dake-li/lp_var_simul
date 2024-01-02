%% DFM SIMULATION STUDY: METHOD CHOICE
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')
addpath(genpath(fullfile('..', 'Subroutines')))
warning('off','MATLAB:structOnObject')

%% SETTINGS

%----------------------------------------------------------------
% Experiment Selection
%----------------------------------------------------------------

% select robustness check mode
mode_select    = 1; % options: 1 (baseline), 2 (small sample), 3 (salient series),
                    % 4 (first diff), 5 (first diff + cumulative IRF),
                    % 6 (first diff + small sample), 7 (first diff + salient series)

% select lag length specifications
lags_select    = 2; % options: 1 (AIC), 2 (4 lags), 3 (8 lags)

% select and group experiments
exper_select_group = {[2,5], [3,6], [1,4]}; % combine G and MP for observed shock, recursive, and IV

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 6 7 8];
methods_obsshock_select  = [1 2 3 4 5 6 7];
methods_recursive_select = [1 2 3 4 5 6 7];

% select a subset of DGPs
select_DGP = 0; % if select a subset of DGPs?
select_DGP_fn = @(i_dgp, res) any(res.settings.specifications.var_select(i_dgp,:)>=res.settings.specifications.random_category_range(11,1)); % binary selection criteria: specifications with asset price & sentiment

% report quantile loss across DGPs in best-method plot
loss_quant = []; % report which quantile loss across DGPs? (default is mean loss, i.e. [])

% Apply shared settings
settings_shared;

%----------------------------------------------------------------
% Preference Plots
%----------------------------------------------------------------

% bias weight grid

n_weight    = 501;
weight_grid = linspace(1,0.5,n_weight)';

% reference method

base_names = {'VAR','LP','BC LP'};
base_indic = NaN(length(exper_files),length(base_names)); % find index of reference method(s)

% construction of choice plots: average over specifications?

choice_averaging = 1;

% lines

lines_plot = lines;
lines_plot = [lines_plot(1:7,:); 0.5 0.5 0.5];
lines_plot = lines_plot([4 3 1 2 8 5 6 7],:);

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

clear cmap

cmap(1,:) = [1 1 1];
cmap(2,:) = [0.5 0.5 0.5];
cmap(3,:) = [0 0 0];

[X,Y] = meshgrid(1:3,1:50);

cmap = interp2(X([1,25,50],:),Y([1,25,50],:),cmap,X,Y);

% turn to coarse grid

n_bin = 5; % Number of bins
length_bin = 50/n_bin;

for i_bin = 1:n_bin
    cmap(1+(i_bin-1)*length_bin:i_bin*length_bin,:) = repmat(cmap(1+(i_bin-1)*length_bin,:),length_bin,1);
end

%% FIGURES

%----------------------------------------------------------------
% Results for Reference Method
%----------------------------------------------------------------

for ne=1:length(exper_files)    
    for nb=1:length(base_names)
        if sum(strcmp(methods_names{ne}, base_names(nb))) == 0
            error('The base method is not included!')
        else
            base_indic(ne,nb) = find(strcmp(methods_names{ne}, base_names(nb)));
        end
    end
end

%----------------------------------------------------------------
% Method Choice Results
%----------------------------------------------------------------

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
            
            % Results relative to true IRF
            the_true_irf = res.DF_model.target_irf; % True IRF
            the_rms_irf = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
            
            the_methods_index = cellfun(@(x) find(strcmp(res.settings.est.methods_name, x)), methods_fields{ne}); % index of each method
            
            the_BIAS2 = extract_struct(res.results.BIAS2);
            the_BIAS2 = the_BIAS2(:,:,the_methods_index);
            the_VCE   = extract_struct(res.results.VCE);
            the_VCE   = the_VCE(:,:,the_methods_index);
    
            the_BIAS2rel = the_BIAS2./the_rms_irf.^2;
            the_VCErel   = the_VCE./the_rms_irf.^2;
            
            %----------------------------------------------------------------
            % Compute Method Choice
            %----------------------------------------------------------------
            
            for nb = 1:length(base_names)
                
                base_method_name = methods_names_plot{base_indic(ne,nb)};
            
                the_objects = methods_names_plot; % Objects to plot
                the_titles =  methods_names_plot; % Plot titles/file names
            
                for j=1:length(the_objects)
                
                    if j == base_indic(ne,nb)
                        continue
                    end
                    
                    % comparison with base method
                    
                    loss_diff = zeros(n_weight,length(res.settings.est.IRF_select),size(the_BIAS2rel,2));
                    for i_weight = 1:n_weight
                        loss_base   = weight_grid(i_weight) * the_BIAS2rel(:,:,base_indic(ne,nb)) + (1-weight_grid(i_weight)) * the_VCErel(:,:,base_indic(ne,nb));
                        loss_method = weight_grid(i_weight) * the_BIAS2rel(:,:,j) + (1-weight_grid(i_weight)) * the_VCErel(:,:,j);
    
                        loss_diff(i_weight,:,:) = loss_method-loss_base;
                    end
                    pref_base = mean(loss_diff>=0,3);
                    
                    if all(abs(loss_diff(:,1,:))<1e-10,'all')
                        the_start_ind=2; % Do not include h=0 if methods are numerically equivalent there
                    else
                        the_start_ind=1;
                    end
                    
                    % plot final results
                    plot_tradeoff(pref_base(:,the_start_ind:end), cmap, horzs(the_start_ind:end)-1, weight_grid, n_bin, font_size);
                    plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_tradeoff_', removeChars(base_method_name), '_vs_', removeChars(the_titles{j}))), output_suffix);
    
                end
            
            end
            
            %----------------------------------------------------------------
            % Compute & Plot Best Overall Procedure
            %----------------------------------------------------------------
            
            if isempty(loss_quant)
                remark_loss_quant = ''; % remark in file name for quantile loss
            else
                remark_loss_quant = strcat('_p', num2str(round(loss_quant*100)));
            end

            choice_raw = zeros(n_weight,max(res.settings.est.IRF_select));
            for i_weight = 1:n_weight

                loss_all = weight_grid(i_weight) * the_BIAS2rel + (1-weight_grid(i_weight)) * the_VCErel;

                if isempty(loss_quant) % mean loss across DGPs
                    loss_all = squeeze(mean(loss_all,2));
                else % quantile loss across DGPs
                    loss_all = squeeze(quantile(loss_all, loss_quant, 2));
                end

                loss_all(:,2:end) = loss_all(:,2:end) + sqrt(eps);
                [~,choice_raw(i_weight,:)] = min(loss_all,[],2);
            end
            
            writematrix([0 horzs; weight_grid(:) choice_raw], fullfile(output_folder, strcat('tradeoff_best', remark_loss_quant, '.csv'))); % Save to file
            
            plot_choice(choice_raw, lines_plot, horzs, weight_grid, methods_select{ne}, ...
                    strjoin({exper_plotname, ': Best Procedure'}), methods_names_plot, 1, font_size);
            plot_save(fullfile(output_folder, strcat(exper_names{ne}, '_tradeoff_best', remark_loss_quant)), output_suffix);  
            
        end
        
    end

end