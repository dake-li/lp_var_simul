%% LP vs VAR: DSGE SIMULATION STUDY
% Dake Li and Christian Wolf
% this version: 3/8/2020

%% HOUSEKEEPING

clc
clear all
close all

path = 'D:\Dake\Princeton\Research\PlagborgMoller_LPVAR\MATLAB_file\Codes';
addpath(genpath([path '/Auxiliary_Functions']))
addpath(genpath([path '/Estimation_Routines']))
cd([path '/DSGE/SW_G']);
addpath(genpath('../Subroutines'))
addpath('C:\Softwares\dynare_4.5.7\matlab')

rng(1);

%% DGP

%----------------------------------------------------------------
% Set up DGP
%----------------------------------------------------------------

% run model

dynare SW_Model noclearall
clean_folder_SW
SW_model.decision = decision(2:end,:);
SW_model.obs = [4 19 5]; % (y,pi,r)

% IV

SW_model.IV.rho     = 0.1; % IV persistence
SW_model.IV.alpha   = 1; % IV shock coefficient
SW_model.IV.sigma_v = 1; % IV noise

%----------------------------------------------------------------
% Represent as ABCDEF Form
%----------------------------------------------------------------

SW_model.n_y   = size(SW_model.obs,2);
SW_model.n_eps = M_.exo_nbr; % 7 structural shocks + instrument noise
SW_model.n_s   = M_.nspred;
SW_model.n_w   = SW_model.n_y;
SW_model.n_e   = SW_model.n_y;

% ABCD representations

SW_model.ABCD = ABCD_fun_SW(SW_model);

% delete superfluous variables

clean_workspace_SW

%% SETTINGS

%----------------------------------------------------------------
% Experiment Specification
%----------------------------------------------------------------

% variable selection

settings.specifications.manual_var_select = 1:1:SW_model.n_y; % select all observables
settings.specifications.random_select     = 0; % no randomization in small model
settings.specifications.plot_indx         = 1; % plot the only specification

% shock position

settings.est.manual_shock_pos        = 1; % manually choose which shock to be our true structural shock?
settings.est.estimate_shock_weight   = 1; % estimate shock weights for true shock? 
settings.est.shock_optimize          = 1; % use shock weights to optimize targeted IRF? otherwise use calibrated shock weights
settings.est.shock_optimize_var      = 1; % for which variable in full model to choose optimal linear combination of shocks 
settings.est.IV_est_normalize_var_pos = 3; % choose IRF normalization variable in all IV methods

% IRFs of interest

settings.est.IRF_response_var_pos = 1; % interested in IRF of which variable in each specification?
settings.est.IRF_hor              = 20; % maximal horizon (include contemporary)
settings.est.IRF_select           = 1:20; % which IRFs to summarize

% compute largest root and VAR(p) fit in population using truncated infinite-order VAR 

settings.est.VAR_infinity_truncate = 50; 

% number of Monte Carlo draws

settings.simul.n_MC    = 10; % number of Monte Carlo reps
settings.simul.seed    = (1:settings.simul.n_MC)*10 + randi([0,9],1,settings.simul.n_MC); % random seed

% variable selections

settings.simul.T      = 200; % time periods for each simulation
settings.simul.T_burn = 100; % burn-in

%----------------------------------------------------------------
% Estimation Settings
%----------------------------------------------------------------

% choose estimand

settings.est.methods_name = {'svar','svar_corrbias','bvar','lp','lp_penalize','var_avg','svar_iv','lp_iv'}; % choose estimands
settings.est.with_shock   = 0; % shock is observed and ordered first in data?
settings.est.recursive_shock = 0; % not use recursive shock
settings.est.with_IV      = 1; % IV is ordered first in data or used in IV method?

% lag specification

settings.est.est_n_lag      = 0; % estimate number of lags?
settings.est.est_n_lag_BIC  = 0; % use BIC? otherwise use AIC
settings.est.n_lags_fix     = 4; % default number of lags if not estimated
settings.est.n_lags_max     = 20; % maximal lag length for info criteria

% BVAR prior

settings.est.prior.tightMN  = 0.1;
settings.est.prior.decay    = 0.5;
settings.est.prior.sig      = 1;
settings.est.prior.tightUR  = 5;
settings.est.prior.tightC   = 5;
settings.est.prior.tightVar = 0.1;

% LP smoothing

settings.est.lambda        = 10; % for lambda = 0 would just do OLS
settings.est.lambdaRange   = [0.005:0.005:0.05, 0.1:0.05:0.5, 1:1:10, 20:10:100]; % cross validation grid, scaled up by T
settings.est.irfLimitOrder = 2; % shrink towards polynomial of that order

% VAR model averaging

settings.est.average_store_weight = [2, 11, 20]; % store model weights at which horizon
settings.est.average_max_lags = 1; % include lags up to n_lags_max? otherwise up to estimated lags

%% PREPARATION

%----------------------------------------------------------------
% Select Specifications
%----------------------------------------------------------------

settings.specifications = pick_var_fn(SW_model, settings);

%----------------------------------------------------------------
% Results Placeholder
%----------------------------------------------------------------

settings.est.n_methods = length(settings.est.methods_name);
settings.est.full_methods_name = {'svar','svar_corrbias','bvar','lp','lp_penalize','var_avg','svar_iv','lp_iv'};

for i_method = 1:length(settings.est.full_methods_name)
    thisMethod = settings.est.methods_name{i_method};
    eval(['results_irf_' thisMethod ...
        '= NaN(settings.est.IRF_hor,settings.simul.n_MC,settings.specifications.n_spec);']); % IRF_hor*n_MC*n_spec
    eval(['results_n_lags_' thisMethod ...
        '= NaN(settings.simul.n_MC,settings.specifications.n_spec);']); %n_MC*n_spec
end
clear i_method thisMethod

results_largest_root_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % n_MC*n_spec
results_LM_stat_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % n_MC*n_spec
results_lambda_lp_penalize = NaN(settings.simul.n_MC,settings.specifications.n_spec); % n_MC*n_spec
results_weight_var_avg = NaN(2*settings.est.n_lags_max,length(settings.est.average_store_weight),...
    settings.simul.n_MC,settings.specifications.n_spec); % n_models*n_horizon*n_MC*n_spec
results_F_stat_svar_iv = NaN(settings.simul.n_MC,settings.specifications.n_spec); %n_MC*n_spec
results_F_stat_lp_iv = NaN(settings.simul.n_MC,settings.specifications.n_spec); %n_MC*n_spec

%% PRELIMINARY COMPUTATIONS: ESTIMANDS

%----------------------------------------------------------------
% Compute True IRFs in Complete Model
%----------------------------------------------------------------

[SW_model.irf, settings.est.shock_weight] = compute_irfs(SW_model,settings);

%----------------------------------------------------------------
% Compute Degree of Invertibility in Specifications
%----------------------------------------------------------------

SW_model.R0_sq = compute_invert_DSGE(SW_model,settings);

%----------------------------------------------------------------
% Compute Persistency of Observables in Specifications
%----------------------------------------------------------------

[SW_model.LRV_Cov_tr_ratio, SW_model.VAR_largest_root, SW_model.frac_coef_for_large_lags] =...
    compute_persist_DSGE(SW_model,settings);

%----------------------------------------------------------------
% Compute True IRF for IV Estimands
%----------------------------------------------------------------

SW_model.IV_irf = compute_IVirfs(SW_model,settings);

%----------------------------------------------------------------
% Compute Popolation IV Strengths
%----------------------------------------------------------------

SW_model.IV_strength = compute_IVstrength_DSGE(SW_model, settings);

%----------------------------------------------------------------
% Compute Target IRF
%----------------------------------------------------------------

SW_model.target_irf = SW_model.IV_irf(settings.est.IRF_select, :);

%% MONTE CARLO ANALYSIS

for i_MC = 1:settings.simul.n_MC

    if mod(i_MC, 10) == 0
        disp("Monte Carlo:")
        disp(i_MC)
    end

    %----------------------------------------------------------------
    % Generate Data
    %----------------------------------------------------------------

    rng(settings.simul.seed(i_MC));

    data_sim_all = generate_data(SW_model,settings);

    %----------------------------------------------------------------
    % List All Temporary Storage for i_MC in parfor
    %----------------------------------------------------------------
    
    temp_irf_svar = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    temp_irf_svar_corrbias = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    temp_irf_bvar = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    temp_irf_lp = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    temp_irf_lp_penalize = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    temp_irf_var_avg = NaN(settings.est.IRF_hor,settings.specifications.n_spec);    
    temp_irf_svar_iv = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    temp_irf_lp_iv = NaN(settings.est.IRF_hor,settings.specifications.n_spec);
    
    temp_n_lags_svar = NaN(1,settings.specifications.n_spec);
    temp_n_lags_svar_corrbias = NaN(1,settings.specifications.n_spec);
    temp_n_lags_bvar = NaN(1,settings.specifications.n_spec);
    temp_n_lags_lp = NaN(1,settings.specifications.n_spec);
    temp_n_lags_lp_penalize = NaN(1,settings.specifications.n_spec);
    temp_n_lags_var_avg = NaN(1,settings.specifications.n_spec);
    temp_n_lags_svar_iv = NaN(1,settings.specifications.n_spec);
    temp_n_lags_lp_iv = NaN(1,settings.specifications.n_spec);
    
    temp_largest_root_svar = NaN(1,settings.specifications.n_spec);
    temp_LM_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_lambda_lp_penalize = NaN(1,settings.specifications.n_spec);
    temp_weight_var_avg = NaN(2*settings.est.n_lags_max,...
        length(settings.est.average_store_weight),settings.specifications.n_spec);
    temp_F_stat_svar_iv = NaN(1,settings.specifications.n_spec);
    temp_F_stat_lp_iv = NaN(1,settings.specifications.n_spec);
    
    %----------------------------------------------------------------
    % Selecting Data
    %----------------------------------------------------------------

    for i_spec = 1:settings.specifications.n_spec
        
        data_sim_select = select_data_fn(data_sim_all,settings,i_spec);
    
        %----------------------------------------------------------------
        % IRF Estimation
        %----------------------------------------------------------------

        % VAR with IV ordered first
        
        if any(strcmp(settings.est.methods_name, 'svar'))
            [temp_irf_svar(:,i_spec),temp_n_lags_svar(1,i_spec),temp_largest_root_svar(1,i_spec),temp_LM_stat_svar(1,i_spec)]...
                = SVAR_est(data_sim_select,settings);
        end

        % bias-corrected VAR with IV ordered first
        
        if any(strcmp(settings.est.methods_name, 'svar_corrbias'))
            [temp_irf_svar_corrbias(:,i_spec),temp_n_lags_svar_corrbias(1,i_spec)]...
                = SVAR_corr_est(data_sim_select,settings);
        end

        % Bayesian VAR with IV ordered first
        
        if any(strcmp(settings.est.methods_name, 'bvar'))
            [temp_irf_bvar(:,i_spec),temp_n_lags_bvar(1,i_spec)]...
                = BVAR_est(data_sim_select,settings);
        end

        % LP with IV ordered first

        if any(strcmp(settings.est.methods_name, 'lp'))
            [temp_irf_lp(:,i_spec),temp_n_lags_lp(1,i_spec)]...
                = LP_est(data_sim_select,settings);
        end

        % shrinkage LP IV ordered first

        if any(strcmp(settings.est.methods_name, 'lp_penalize'))
            [temp_irf_lp_penalize(:,i_spec),temp_n_lags_lp_penalize(1,i_spec), temp_lambda_lp_penalize(1,i_spec)]...
                = LP_shrink_est(data_sim_select,settings);
        end

        % VAR model averaging with IV ordered first
        
        if any(strcmp(settings.est.methods_name, 'var_avg'))
            [temp_irf_var_avg(:,i_spec),temp_n_lags_var_avg(1,i_spec), temp_weight_var_avg(:,:,i_spec)]...
                = VAR_avg_est(data_sim_select,settings);
        end
        
        % SVAR-IV       

        if any(strcmp(settings.est.methods_name, 'svar_iv'))
            [temp_irf_svar_iv(:,i_spec),temp_n_lags_svar_iv(1,i_spec),temp_F_stat_svar_iv(1,i_spec)]...
                = SVAR_IV_est(data_sim_select,settings);
        end

        % LP-IV

        if any(strcmp(settings.est.methods_name, 'lp_iv'))
            [temp_irf_lp_iv(:,i_spec),temp_n_lags_lp_iv(1,i_spec),temp_F_stat_lp_iv(1,i_spec)]...
                = LP_IV_est(data_sim_select,settings);
        end
            
    end
    
    %----------------------------------------------------------------
    % Move Results to Permanent Storage in parfor
    %----------------------------------------------------------------
    
    results_irf_svar(:,i_MC,:) = temp_irf_svar;
    results_irf_svar_corrbias(:,i_MC,:) = temp_irf_svar_corrbias;
    results_irf_bvar(:,i_MC,:) = temp_irf_bvar;
    results_irf_lp(:,i_MC,:) = temp_irf_lp;
    results_irf_lp_penalize(:,i_MC,:) = temp_irf_lp_penalize;
    results_irf_var_avg(:,i_MC,:) = temp_irf_var_avg;
    results_irf_svar_iv(:,i_MC,:) = temp_irf_svar_iv;
    results_irf_lp_iv(:,i_MC,:) = temp_irf_lp_iv;
    
    results_n_lags_svar(i_MC,:) = temp_n_lags_svar;
    results_n_lags_svar_corrbias(i_MC,:) = temp_n_lags_svar_corrbias;
    results_n_lags_bvar(i_MC,:) = temp_n_lags_bvar;
    results_n_lags_lp(i_MC,:) = temp_n_lags_lp;
    results_n_lags_lp_penalize(i_MC,:) = temp_n_lags_lp_penalize;
    results_n_lags_var_avg(i_MC,:) = temp_n_lags_var_avg;
    results_n_lags_svar_iv(i_MC,:) = temp_n_lags_svar_iv;
    results_n_lags_lp_iv(i_MC,:) = temp_n_lags_lp_iv;
    
    results_largest_root_svar(i_MC,:) = temp_largest_root_svar;
    results_LM_stat_svar(i_MC,:) = temp_LM_stat_svar;
    results_lambda_lp_penalize(i_MC,:) = temp_lambda_lp_penalize;
    results_weight_var_avg(:,:,i_MC,:) = temp_weight_var_avg;
    results_F_stat_svar_iv(i_MC,:) = temp_F_stat_svar_iv;
    results_F_stat_lp_iv(i_MC,:) = temp_F_stat_lp_iv;

end

% clear temporary storage

for i_method = 1:length(settings.est.full_methods_name)
    thisMethod = settings.est.full_methods_name{i_method};
    eval(['clear temp_irf_' thisMethod ';']);
    eval(['clear temp_n_lags_' thisMethod ';']);
end
clear temp_largest_root_svar temp_LM_stat_svar temp_lambda_lp_penalize temp_weight_var_avg temp_F_stat_svar_iv temp_F_stat_lp_iv
clear i_MC i_spec data_sim_all data_sim_select i_method thisMethod

%% SUMMARIZE RESULTS

%----------------------------------------------------------------
% Wrap up Results, Pick out Target IRF
%----------------------------------------------------------------

% wrap up results from parallel loop

for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    eval(['results.irf.' thisMethod '= results_irf_' thisMethod ';']);
    eval(['results.n_lags.' thisMethod '= results_n_lags_' thisMethod ';']);
    
end

if any(strcmp(settings.est.methods_name, 'svar'))    
    results.largest_root.svar = results_largest_root_svar;
    results.LM_stat.svar = results_LM_stat_svar;
end

if any(strcmp(settings.est.methods_name, 'var_avg'))
    results.weight.var_avg = results_weight_var_avg;
end

if any(strcmp(settings.est.methods_name, 'svar_iv'))
    results.F_stat.svar_iv = results_F_stat_svar_iv;
end

if any(strcmp(settings.est.methods_name, 'lp_iv'))
    results.F_stat.lp_iv = results_F_stat_lp_iv;
end

for i_method = 1:length(settings.est.full_methods_name)
    thisMethod = settings.est.full_methods_name{i_method};
    eval(['clear results_irf_' thisMethod ' results_n_lags_' thisMethod ';']);
end

clear results_largest_root_svar results_LM_stat_svar results_lambda_lp_penalize results_weight_var_avg results_F_stat_svar_iv results_F_stat_lp_iv
clear i_method thisMethod

% store IRF only at a few horizons

for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    eval(['results.irf.' thisMethod '= results.irf.' thisMethod '(settings.est.IRF_select,:,:);']);
    
end

clear i_method thisMethod

%----------------------------------------------------------------
% Compute Mean-Squared Errors, Bias-Squared, Variance
%----------------------------------------------------------------

% compute MSE, Bias2, VCE for each horizon and each specification

for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    
    eval(['results.MSE.' thisMethod '= squeeze(mean((results.irf.' thisMethod ...
        ' - permute(SW_model.target_irf,[1 3 2])).^2, 2));']);
    
    eval(['results.BIAS2.' thisMethod '= (squeeze(mean(results.irf.' thisMethod ...
        ', 2)) - SW_model.target_irf).^2;']);
    
    eval(['results.VCE.' thisMethod '= squeeze(var(results.irf.' thisMethod ', 0, 2));']);
    
end

clear i_method thisMethod

%% PLOT RESULTS

%----------------------------------------------------------------
% Plot IRFs for Checking
%----------------------------------------------------------------

% for i_method = 1:settings.est.n_methods
%     
%     thisMethod = settings.est.methods_name{i_method};
%     figure(i_method)
%     plot(settings.est.IRF_select, SW_model.target_irf(:,settings.specifications.plot_indx),'Linewidth',5)
%     hold on
%     for i = 1:settings.simul.n_MC
%         eval(['plot(settings.est.IRF_select, results.irf.' thisMethod '(:,i,settings.specifications.plot_indx))']);
%         hold on
%     end
%     title(replace(thisMethod,'_',' '))
%     hold off
% 
% end
% 
% clear i_method thisMethod i