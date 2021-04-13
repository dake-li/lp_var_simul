%% DFM SIMULATION STUDY: PLOT BIAS/VARIANCE FRONTIER
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

%% HOUSEKEEPING

clc
clear all
close all

addpath('Plotting_Functions')

%% SETTINGS

% select lag length specifications
lags_select    = 2; % options: 1 (AIC), 2 (4 lags), 3 (8 lags)

% select and group experiments
exper_select_group = {[2,5]}; % combine G and MP for observed shock, recursive, and IV

% select estimation methods for each experiment
methods_iv_select        = [1 2 3 4 5 6 7];
methods_obsshock_select  = [1 2 3 4 5 6];
methods_recursive_select = [1 2 3 4 5 6];

% Apply shared settings
settings_shared;

%% FIGURES

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
        
        %----------------------------------------------------------------
        % Compute Reporting Results
        %----------------------------------------------------------------
        
        the_true_irf = res.DF_model.target_irf; % True IRF
        the_rms_irf  = sqrt(mean(the_true_irf.^2)); % Root average squared true IRF across horizons
        
        % frontier settings
        
        frontier_hor_lb = 4;
        frontier_hor_ub = 16;
        
        % bias
        
        the_bias = squeeze(median(sqrt(extract_struct(res.results.BIAS2))./the_rms_irf, 2));
        the_bias = mean(the_bias(frontier_hor_lb:frontier_hor_ub,:),1)';
        
        % standard deviation
        
        the_std = squeeze(median(sqrt(extract_struct(res.results.VCE))./the_rms_irf, 2));        
        the_std = mean(the_std(frontier_hor_lb:frontier_hor_ub,:),1)';
        
        % approximate frontier
        
        bias_frontier = linspace(0.25 * min(the_bias),1.25 * max(the_bias),20)';
        
        obj_fn    = @(poly) obj_fn_aux(poly,the_bias,the_std);
        constr_fn = @(poly) constr_fn_aux(poly,the_bias,the_std,bias_frontier);
        
        solution = fmincon(obj_fn,[0 0 0],[],[],[],[],[],[],constr_fn);
%         solution = fmincon(obj_fn,[0 0 0],[],[],[],[],[],[],[]);
        
        std_frontier  = solution(1) + solution(2) * bias_frontier + solution(3) * bias_frontier.^2;
        
        %----------------------------------------------------------------
        % Plot Results
        %----------------------------------------------------------------
        
        plot_frontier([the_bias,the_std],[bias_frontier,std_frontier],methods_names_plot,20);
        
    end
    
end