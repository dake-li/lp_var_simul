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
        
        frontier_hor_lb_1 = 3;
        frontier_hor_ub_1 = 3;
        frontier_hor_lb_2 = 8;
        frontier_hor_ub_2 = 8;
        frontier_hor_lb_3 = 4;
        frontier_hor_ub_3 = 16;
        
        % bias
        
        the_bias = squeeze(median(sqrt(extract_struct(res.results.BIAS2))./the_rms_irf, 2));
        
        the_bias_1 = mean(the_bias(frontier_hor_lb_1:frontier_hor_ub_1,:),1)';
        the_bias_2 = mean(the_bias(frontier_hor_lb_2:frontier_hor_ub_2,:),1)';
        the_bias_3 = mean(the_bias(frontier_hor_lb_3:frontier_hor_ub_3,:),1)';
        
        % standard deviation
        
        the_std = squeeze(median(sqrt(extract_struct(res.results.VCE))./the_rms_irf, 2));  
        
        the_std_1 = mean(the_std(frontier_hor_lb_1:frontier_hor_ub_1,:),1)';
        the_std_2 = mean(the_std(frontier_hor_lb_2:frontier_hor_ub_2,:),1)';
        the_std_3 = mean(the_std(frontier_hor_lb_3:frontier_hor_ub_3,:),1)';
        
        % approximate frontier
        
        bias_frontier_1 = linspace(0.25 * min(the_bias_1),1.25 * max(the_bias_1),20)';
        bias_frontier_2 = linspace(0.25 * min(the_bias_2),1.25 * max(the_bias_2),20)';
        bias_frontier_3 = linspace(0.25 * min(the_bias_3),1.25 * max(the_bias_3),20)';
        
        std_frontier_1 = std_frontier_fn(the_bias_1,the_std_1,bias_frontier_1);
        std_frontier_2 = std_frontier_fn(the_bias_2,the_std_2,bias_frontier_2);
        std_frontier_3 = std_frontier_fn(the_bias_3,the_std_3,bias_frontier_3);
        
        %----------------------------------------------------------------
        % Plot Results
        %----------------------------------------------------------------
        
        pos_1  = [2 3 1 1 1 1];
        dist_1 = [0.05 1 0.25 0.25 0.25 0.25];       
        plot_frontier([the_bias_1,the_std_1],[bias_frontier_1,std_frontier_1],methods_names_plot,20,...
            ['Horizon: h = ' num2str(frontier_hor_lb_1)],pos_1,dist_1);
        
        pos_2  = [1 3 1 1 1 1];
        dist_2 = [0.05 1 0.25 0.25 0.25 0.25];       
        plot_frontier([the_bias_2,the_std_2],[bias_frontier_2,std_frontier_2],methods_names_plot,20,...
            ['Horizon: h = ' num2str(frontier_hor_lb_2)],pos_2,dist_2);
        
        pos_3  = [1 3 1 1 1 1];
        dist_3 = [0.05 1 0.25 0.25 0.25 0.25];       
        plot_frontier([the_bias_3,the_std_3],[bias_frontier_3,std_frontier_3],methods_names_plot,20,...
            ['Average between h = ' num2str(frontier_hor_lb_3) ' and h = ' num2str(frontier_hor_ub_3)],pos_3,dist_3);
        
    end
    
end