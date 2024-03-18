%% DFM SIMULATION STUDY: MAIN FILE
% Dake Li, Mikkel Plagborg-M�ller and Christian Wolf

%% HOUSEKEEPING

clc
clear all
close all

addpath(genpath(fullfile('..', 'Auxiliary_Functions')))
addpath(genpath(fullfile('..', 'Estimation_Routines')))
addpath(genpath('Subroutines'))

rng(1, 'twister');
tic;

% Parallel computing object

num_workers = 1; % number of workers in a parallel pool
poolobj = parpool('local', num_workers);
clear num_workers;

%% SET EXPERIMENT

spec_id = 1; % seed for random draws of specifications (= DGPs from encompassing model)
dgp_type = 'G'; % structural shock: either 'G' or 'MP'
estimand_type = 'ObsShock'; % structural estimand: either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (= AIC)
mode_type = 1; % robustness check mode:
               % 1 (baseline), 2 (small sample), 3 (large sample),
               % 4 (salient series), 5 (more observables), 6 (first diff)

estim_diagn = 0; % =1: show DFM estimation diagnostics

%% SETTINGS

% Apply shared settings as well as settings specific to DGP and estimand type

run(fullfile('Settings', 'shared'));
run(fullfile('Settings', dgp_type));
run(fullfile('Settings', estimand_type));
run(fullfile('Settings', 'check_mode'));

% Storage folder for results

save_pre = 'Results'; % destination to store the results

if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));

%% ENCOMPASSING DFM MODEL

%----------------------------------------------------------------
% Set up Encompassing Model
%----------------------------------------------------------------

% estimate DFM from dataset

DFM_estimate = DFM_est(DF_model.n_fac, DF_model.n_lags_fac, DF_model.n_lags_uar, DF_model.reorder, DF_model.levels, DF_model.coint_rank);
if estim_diagn == 1
    run_estim_diagn; % Display several estimation diagnostics
end

% extract and store estimated DFM parameters

DF_model.Phi           = DFM_estimate.Phi;
DF_model.Sigma_eta     = DFM_estimate.Sigma_eta;

DF_model.Lambda        = DFM_estimate.Lambda;
DF_model.delta         = DFM_estimate.delta;
DF_model.sigma_v       = DFM_estimate.sigma_v;

DF_model.variable_name_code = DFM_estimate.bpnamevec;
DF_model.variable_name_short = DFM_estimate.bplabvec_short;
DF_model.variable_name_long = DFM_estimate.bplabvec_long;
DF_model.trans_code = DFM_estimate.bptcodevec; % transformation code
% (1) y = x, (2) y = (1-L)x, (3) y = (1-L)^2 x,
% (4) y = ln(x), (5) y = (1-L)ln(x), (6) y = (1-L)^2 ln(x)

%----------------------------------------------------------------
% Set Up IV DGP
%----------------------------------------------------------------

if strcmp(estimand_type, 'IV')

    DF_model.IV.rho_grid = DF_model.IV.rho * settings.est.IV.IV_persistence_scale;
    DF_model.IV.sigma_v_grid = DF_model.IV.sigma_v * settings.est.IV.IV_strength_scale;

end

%----------------------------------------------------------------
% Represent as Model in ABCD Form
%----------------------------------------------------------------

[DF_model.n_y,DF_model.n_fac] = size(DF_model.Lambda);
DF_model.ABCD  = ABCD_fun_DFM(DF_model);

%----------------------------------------------------------------
% Shock Weights
%----------------------------------------------------------------

shock_weight = zeros(DF_model.n_fac+DF_model.n_y,1);
if settings.est.estimate_shock_weight==1 % if want to estimate the optimal shock weight to maximize impact response
    shock_weight(1:DF_model.n_fac) = DF_model.ABCD.D(settings.est.shock_optimize_var_IRF,1:DF_model.n_fac); % optimal weight is in the same direction as impulse response
    shock_weight = shock_weight/sqrt(shock_weight'*shock_weight); % normalize weight
else
    shock_weight = zeros(n_eps, 1); % manually choose shock
    shock_weight(settings.est.manual_shock_pos) = 1;
end
settings.est.shock_weight = shock_weight;
clear shock_weight;

%% PREPARATION

%----------------------------------------------------------------
% Select Individual DGPs from Encompassing Model
%----------------------------------------------------------------

% randomly draw DGPs

settings.specifications = pick_var_fn(DF_model, settings, spec_id);

%----------------------------------------------------------------
% Create Placeholders for Results
%----------------------------------------------------------------

% number of estimation methods

settings.est.n_methods = length(settings.est.methods_name);

% impulse response estimates

results_irf = NaN(settings.est.n_methods,settings.est.IRF_hor,settings.simul.n_MC,settings.specifications.n_spec); % estimated IRFs for each method: size IRF_hor*n_MC*n_spec

% several other features of each Monte Carlo estimation

results_n_lags = NaN(settings.est.n_methods,settings.simul.n_MC,settings.specifications.n_spec); % estimated lags for each method: size n_MC*n_spec
results_largest_root_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % largest VAR root: size n_MC*n_spec
results_LM_stat_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % LM statistic: size n_MC*n_spec
results_LM_pvalue_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % LM p value: size n_MC*n_spec
results_Hausman_stat_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % LM statistic: size n_MC*n_spec
results_Hausman_pvalue_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % LM p value: size n_MC*n_spec
results_Granger_stat_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % Granger statistic: size n_MC*n_spec
results_Granger_pvalue_svar = NaN(settings.simul.n_MC,settings.specifications.n_spec); % Granger p value: size n_MC*n_spec
results_GLP_hyper_bvar = NaN(3,settings.simul.n_MC,settings.specifications.n_spec); % GLP hyperparameters: size 3*n_MC*n_spec
results_lambda_lp_penalize = NaN(settings.simul.n_MC,settings.specifications.n_spec); % pen. LP lambda: size n_MC*n_spec
results_weight_var_avg = NaN(2*settings.est.n_lags_max,length(settings.est.average_store_weight),...
    settings.simul.n_MC,settings.specifications.n_spec); % weights in VAR averaging: size n_models*n_horizon*n_MC*n_spec
results_submodel_irf_var_avg = NaN(2*settings.est.n_lags_max,settings.est.IRF_hor,...
    settings.simul.n_MC,settings.specifications.n_spec); % IRF in each VAR submodel: size n_models*IRF_hor*n_MC*n_spec
results_F_stat_svar_iv = NaN(settings.simul.n_MC,settings.specifications.n_spec); % VAR-IV F statistic: size n_MC*n_spec
results_F_pvalue_svar_iv = NaN(settings.simul.n_MC,settings.specifications.n_spec); % VAR-IV F p value: size n_MC*n_spec

%% STRUCTURAL ESTIMANDS AND DGP SUMMARIES

DF_model = dgp_irfs_stats(DF_model, settings, estimand_type);

%% MONTE CARLO ANALYSIS

disp('Monte Carlo simulation starts.');
disp(['specification choice set id: ', num2str(spec_id)]);
disp(['dgp type: ', dgp_type]);
disp(['estimand type: ', estimand_type]);
disp(['lag type: ', num2str(lag_type)]);

parfor i_MC = 1:settings.simul.n_MC

    if mod(i_MC, 100) == 0
        disp(['Monte Carlo repetitions: ', num2str(i_MC)])
    end

    %----------------------------------------------------------------
    % Generate Data From Encompassing DFM
    %----------------------------------------------------------------
    
    rng(settings.simul.seed(i_MC), 'twister');

    data_sim_all = generate_data(DF_model,settings);

    %----------------------------------------------------------------
    % Several Temporary Storage Folders (due to parfor)
    %----------------------------------------------------------------
    
    temp_irf = NaN(settings.est.n_methods,settings.est.IRF_hor,settings.specifications.n_spec);
    temp_n_lags = NaN(settings.est.n_methods,settings.specifications.n_spec);
    
    temp_largest_root_svar = NaN(1,settings.specifications.n_spec);
    temp_LM_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_LM_pvalue_svar = NaN(1,settings.specifications.n_spec);
    temp_Hausman_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_Hausman_pvalue_svar = NaN(1,settings.specifications.n_spec);
    temp_Granger_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_Granger_pvalue_svar = NaN(1,settings.specifications.n_spec);
    temp_GLP_hyper_bvar = NaN(3,1,settings.specifications.n_spec);
    temp_lambda_lp_penalize = NaN(1,settings.specifications.n_spec);
    temp_weight_var_avg = NaN(2*settings.est.n_lags_max,...
        length(settings.est.average_store_weight),settings.specifications.n_spec);
    temp_submodel_irf_var_avg = NaN(2*settings.est.n_lags_max,...
        settings.est.IRF_hor,settings.specifications.n_spec);
    temp_F_stat_svar_iv = NaN(1,settings.specifications.n_spec);
    temp_F_pvalue_svar_iv = NaN(1,settings.specifications.n_spec);
    
    %----------------------------------------------------------------
    % Estimation for each DGP
    %----------------------------------------------------------------

    for i_spec = 1:settings.specifications.n_spec
        
        % select data
        
        data_sim_select = select_data_fn(data_sim_all,settings,i_spec);
    
        % estimate IRFs
        
        for i_method = 1:settings.est.n_methods
            
            switch settings.est.methods_name{i_method}

                case 'svar' % VAR
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec),temp_largest_root_svar(i_spec),...
                        temp_LM_stat_svar(i_spec),temp_LM_pvalue_svar(i_spec),...
                        temp_Hausman_stat_svar(i_spec),temp_Hausman_pvalue_svar(i_spec),...
                        temp_Granger_stat_svar(i_spec),temp_Granger_pvalue_svar(i_spec)]...
                        = SVAR_est(data_sim_select,settings,0);

                case 'svar_corrbias' % bias-corrected VAR
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec)]...
                        = SVAR_est(data_sim_select,settings,1);

                case 'bvar' % Bayesian VAR
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec),temp_GLP_hyper_bvar(:,1,i_spec)]...
                        = BVAR_est(data_sim_select,settings);

                case 'lp' % LP
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec)]...
                        = LP_est(data_sim_select,settings,0);

                case 'lp_corrbias' % bias-corrected LP
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec)]...
                        = LP_est(data_sim_select,settings,1);

                case 'lp_penalize' % shrinkage LP
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec), temp_lambda_lp_penalize(i_spec)]...
                        = LP_shrink_est(data_sim_select,settings);

                case 'var_avg' % VAR model averaging
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec), temp_weight_var_avg(:,:,i_spec), temp_submodel_irf_var_avg(:,:,i_spec)]...
                        = VAR_avg_est(data_sim_select,settings);

                case 'svar_iv' % SVAR-IV       
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec),...
                        temp_F_stat_svar_iv(i_spec),temp_F_pvalue_svar_iv(i_spec)]...
                        = SVAR_IV_est(data_sim_select,settings);
                
            end
            
        end
        
    end
    
    %----------------------------------------------------------------
    % Move Results to Permanent Storage in parfor
    %----------------------------------------------------------------
    
    results_irf(:,:,i_MC,:) = temp_irf;
    results_n_lags(:,i_MC,:) = temp_n_lags;
    
    results_largest_root_svar(i_MC,:) = temp_largest_root_svar;
    results_LM_stat_svar(i_MC,:) = temp_LM_stat_svar;
    results_LM_pvalue_svar(i_MC,:) = temp_LM_pvalue_svar;
    results_Hausman_stat_svar(i_MC,:) = temp_Hausman_stat_svar;
    results_Hausman_pvalue_svar(i_MC,:) = temp_Hausman_pvalue_svar;
    results_Granger_stat_svar(i_MC,:) = temp_Granger_stat_svar;
    results_Granger_pvalue_svar(i_MC,:) = temp_Granger_pvalue_svar;
    results_GLP_hyper_bvar(:,i_MC,:) = temp_GLP_hyper_bvar;
    results_lambda_lp_penalize(i_MC,:) = temp_lambda_lp_penalize;
    results_weight_var_avg(:,:,i_MC,:) = temp_weight_var_avg;
    results_submodel_irf_var_avg(:,:,i_MC,:) = temp_submodel_irf_var_avg;
    results_F_stat_svar_iv(i_MC,:) = temp_F_stat_svar_iv;
    results_F_pvalue_svar_iv(i_MC,:) = temp_F_pvalue_svar_iv;

end

% clear temporary storage

clear temp_* i_MC i_spec data_sim_all data_sim_select i_method

%% SUMMARIZE RESULTS

%----------------------------------------------------------------
% Pick Out Main Results
%----------------------------------------------------------------

% extract IRF and lag results for all methods

for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    results.irf.(thisMethod) = permute(results_irf(i_method,settings.est.IRF_select,:,:), [2 3 4 1]);
    results.n_lags.(thisMethod) = permute(results_n_lags(i_method,:,:), [2 3 1]);
    
end

% extract additional method-specific results

if any(strcmp(settings.est.methods_name, 'svar'))    
    results.largest_root.svar = results_largest_root_svar;
    results.LM_stat.svar = results_LM_stat_svar;
    results.LM_pvalue.svar = results_LM_pvalue_svar;
    if strcmp(estimand_type, 'ObsShock')
        results.Hausman_stat.svar = results_Hausman_stat_svar;
        results.Hausman_pvalue.svar = results_Hausman_pvalue_svar;
    end
    if strcmp(estimand_type, 'IV')
        results.Granger_stat.svar = results_Granger_stat_svar;
        results.Granger_pvalue.svar = results_Granger_pvalue_svar;
    end
end

if any(strcmp(settings.est.methods_name, 'bvar'))
    if settings.est.bvar_glp == 1
        results.GLP_hyper.lambda = squeeze(results_GLP_hyper_bvar(1,:,:)); % Std dev of MN prior
        results.GLP_hyper.delta = squeeze(results_GLP_hyper_bvar(2,:,:)); % Std dev of dummy-initial-observation prior
        results.GLP_hyper.mu = squeeze(results_GLP_hyper_bvar(3,:,:)); % Std dev of sum-of-coefficients prior
    end
end

if any(strcmp(settings.est.methods_name, 'lp_penalize'))    
    results.lambda.lp_penalize = results_lambda_lp_penalize;
end

if any(strcmp(settings.est.methods_name, 'var_avg'))
    results.weight.var_avg = results_weight_var_avg;
    if settings.est.average_store_submodel_irf == 1
        results.submodel_irf.var_avg = results_submodel_irf_var_avg;
    end
end

if any(strcmp(settings.est.methods_name, 'svar_iv'))
    results.F_stat.svar_iv = results_F_stat_svar_iv;
    results.F_pvalue.svar_iv = results_F_pvalue_svar_iv;
end

clear results_* i_method thisMethod

%----------------------------------------------------------------
% Export Results
%----------------------------------------------------------------

mkdir(save_folder);
save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type, '_', num2str(spec_id))), ...
    'DFM_estimate','DF_model','settings','results',...
    'spec_id','dgp_type','estimand_type','lag_type', 'mode_type', '-v7.3'); % save results

delete(poolobj);
clear save_folder save_pre save_mode_dir save_suff poolobj mode_list

toc;