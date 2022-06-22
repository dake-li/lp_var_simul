%% PLOT TRUE IRF WITH SAMPLE IRF ESTIMATES
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clear all;
addpath(genpath(fullfile('..', 'Reporting', 'Plotting_Functions')));

%% SET UP DESTINATION FOLDER AND FILES

% Experiment setup

spec_id = 1; % specification choice set id
dgp_type = 'G'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)
mode_type = 1; % robustness check mode:
               % 1 (baseline), 2 (cumulative IRF), 
               % 3 (persistent DGP), 4 (persistent DGP with MN prior), 
               % 5 (small sample), 6 (salient series)

% Select methods

methods_select = 1:6; % choose index from ('VAR','BC VAR','BVAR','LP','Pen LP','VAR Avg','SVAR-IV')

% Plot settings

rng(1, 'twister');
dgp_idx_to_plot = 3; % plot which DGP
n_estimate_randomly_select = 10; % randomly select several estimates to plot

font_size = 10;
output_suffix = 'png';

% Storage folder for results

save_pre = fullfile('..', 'Results');

mode_list   = {'baseline', 'cumulative', 'persistent', 'persistent_BVAR_MN_prior' , 'small', 'salient'};
save_mode_dir = mode_list{mode_type}; % set up directory for robustness-check modes

if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));

%% PLOT IRF

% load results

filename = fullfile(save_folder, strcat(strcat('DFM_', dgp_type, '_', estimand_type), '_', num2str(spec_id)));
res = load(filename);

% plot true and estimated irf for each method

n_estimate_used = min(n_estimate_randomly_select, res.settings.simul.n_MC); % number of estimates finally used
MC_idx_set = randsample(1:res.settings.simul.n_MC, n_estimate_used); % which Monte Carlo estimates to plot

methods_fields       = {'svar','svar_corrbias','bvar','lp','lp_penalize','var_avg','svar_iv'};
methods_names        = {'VAR','BC VAR','BVAR','LP','Pen LP','VAR Avg','SVAR-IV'};

for i_method = methods_select
    
    thisMethodField = methods_fields{i_method};
    thisMethodName = methods_names{i_method};
    
    % setup plot

    figure;
    set(gca,'FontSize', font_size);
    set(gca,'TickLabelInterpreter','latex')
    grid on;
    hold on;

    % true irf

    plot(res.settings.est.IRF_select-1, res.DF_model.target_irf(:, dgp_idx_to_plot),'Linewidth', 5);
    
    
    % estimated irf

    for i_MC = MC_idx_set
        plot(res.settings.est.IRF_select-1, res.results.irf.(thisMethodField)(:, i_MC, dgp_idx_to_plot), '--','LineWidth', 1);
        hold on
    end

    % add label

    hold off;
    % title(strjoin({estimand_type, dgp_type, ':', thisMethodName, 'IRFs'}));
    xlabel('Horizon','interpreter','latex','FontSize',font_size);
    set(gca, 'XTick', [min(res.settings.est.IRF_select-1) 2:2:max(res.settings.est.IRF_select-1)]);

    % save plot

    plot_save(fullfile(save_folder, strcat('irf_estimate_', removeChars(thisMethodName))), output_suffix);
 
end