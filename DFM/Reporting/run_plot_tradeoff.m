%% PLOT SIMULATION RESULTS: BIAS/VARIANCE TRADE-OFF
% this version: 09/22/2020

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')
warning('off','MATLAB:structOnObject')

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

%----------------------------------------------------------------
% Preference Plots
%----------------------------------------------------------------

% bias weight grid

n_weight    = 10001;
weight_grid = linspace(1,0,n_weight)';

% reference method (should have something here on this being LP...)

% base_names = {'SVAR','LP'};
base_names = {'SVAR'};

% bias-variance ordering

trade_off_pos_iv = [4 3 6 1 2 5 7];   % IV DGP
trade_off_pos_noniv = [4 3 6 1 2 5];  % Other DGPs

% construction of choice plots: average over specifications?

choice_averaging = 1;

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

n = 200;
clear cmap

cmap(1,:) = [1 1 1];
cmap(2,:) = [0.5 0.5 0.5];
cmap(3,:) = [0 0 0];

[X,Y] = meshgrid([1:3],[1:50]);

cmap = interp2(X([1,25,50],:),Y([1,25,50],:),cmap,X,Y);

cmap_inv(1,:) = [0 0 0];
cmap_inv(2,:) = [0.5 0.5 0.5];
cmap_inv(3,:) = [1 1 1];

cmap_inv = interp2(X([1,25,50],:),Y([1,25,50],:),cmap_inv,X,Y);

%% FIGURES

% Find index of reference method
base_indic = NaN(length(exper_files),length(base_names));

for ne=1:length(exper_files)
    for nb=1:length(base_names)
        if sum(strcmp(methods_names{ne}, base_names(nb))) == 0
            error('The base method is not included!')
        else
            base_indic(ne,nb) = find(strcmp(methods_names{ne}, base_names(nb)));
        end
    end
end

for nf=1:length(lags_folders) % For each folder...

    for ne=1:length(exper_files) % For each experiment in folder...
        
        % Load results
        load_results;
        
        % Results relative to true IRF
        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        the_BIAS2 = extract_struct(res.results.BIAS2);
        the_BIAS2 = the_BIAS2(:,:,methods_select{ne});
        the_VCE   = extract_struct(res.results.VCE);
        the_VCE   = the_VCE(:,:,methods_select{ne});

        the_BIAS2rel = the_BIAS2./the_rms_irf.^2;
        the_VCErel   = the_VCE./the_rms_irf.^2;
        
        %----------------------------------------------------------------
        % Compute Comparison Results
        %----------------------------------------------------------------
        
        for nb = 1:length(base_names)
        
        base_method_name = methods_names_plot{base_indic(ne,nb)};
        
        the_objects = methods_names_plot; % Objects to plot
        the_titles =  methods_names_plot; % Plot titles/file names

        for j=1:length(the_objects)
            
            if j == base_indic(ne,nb)
                continue
            end
            
            pref_base = zeros(n_weight,max(res.settings.est.IRF_select),size(the_BIAS2rel,2));
            for i_weight = 1:n_weight
                loss_base   = weight_grid(i_weight) * the_BIAS2rel(:,:,base_indic(ne,nb)) + (1-weight_grid(i_weight)) * the_VCErel(:,:,base_indic(ne,nb));
                loss_method = weight_grid(i_weight) * the_BIAS2rel(:,:,j) + (1-weight_grid(i_weight)) * the_VCErel(:,:,j);

                pref_base(i_weight,:,:) = (loss_base <= loss_method);
            end
            pref_base = mean(pref_base,3);

            plot_tradeoff(1-pref_base(:,2:end), cmap, horzs(2:end)-1, weight_grid, ...
                strjoin({exper_plotname, ':', the_titles{j}, 'Preferred Over', base_method_name}), font_size)
            plot_save(fullfile(output_folder, strcat('tradeoff_', removeChars(the_titles{j}), '_vs_', removeChars(base_method_name))), output_suffix);

        end
        
        end
        
        %----------------------------------------------------------------
        % Compute Best Procedure
        %----------------------------------------------------------------
        
        choice_raw = zeros(n_weight,max(res.settings.est.IRF_select));
        for i_weight = 1:n_weight
            loss_all = weight_grid(i_weight) * the_BIAS2rel + (1-weight_grid(i_weight)) * the_VCErel;
            loss_all = squeeze(mean(loss_all,2));
            [~,choice_raw(i_weight,:)] = min(loss_all,[],2);
        end
        
        plot_choice(choice_raw(:,2:end), lines, horzs(2:end)-1, weight_grid, methods_select{ne}, ...
                strjoin({exper_plotname, ': Best Procedure'}), methods_names_plot, 1, font_size);
        plot_save(fullfile(output_folder, 'tradeoff_best'), output_suffix);    
        
        %----------------------------------------------------------------
        % Compute Choice Results
        %----------------------------------------------------------------
        
        if isempty(strfind(exper_names{ne}, 'IV'))
            trade_off_pos = trade_off_pos_noniv;
        else
            trade_off_pos = trade_off_pos_iv;
        end
        [~,trade_off_aux] = sort(trade_off_pos);
        
        if choice_averaging == 0
            
            choice = NaN(n_weight,max(res.settings.est.IRF_select));
            for i_method = 1:length(methods_names{ne})
                choice(choice_raw == i_method) = trade_off_pos(methods_select{ne}(i_method));
            end
        
        else

            choice_raw = zeros(n_weight,max(res.settings.est.IRF_select),size(the_BIAS2rel,2));
            for i_weight = 1:n_weight
                loss_all = weight_grid(i_weight) * the_BIAS2rel + (1-weight_grid(i_weight)) * the_VCErel;
                [~,choice_raw(i_weight,:,:)] = min(loss_all,[],3);
            end
            choice = NaN(n_weight,max(res.settings.est.IRF_select),size(the_BIAS2rel,2));
            for i_method = 1:length(methods_names{ne})
                choice(choice_raw == i_method) = trade_off_pos(methods_select{ne}(i_method));
            end
            choice = mean(choice,3);
        
        end
        
        plot_choice(choice(:,2:end), cmap_inv, horzs(2:end)-1, weight_grid, methods_select{ne}, ...
                strjoin({exper_plotname, ': Method Choice'}), methods_iv_names(trade_off_aux), 0, font_size);
        plot_save(fullfile(output_folder, 'tradeoff_choice'), output_suffix);
        
    end
    
end