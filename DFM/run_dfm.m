%% LP vs VAR: DFM SIMULATION STUDY
% Dake Li and Christian Wolf


%% HOUSEKEEPING

clc
clear all
close all

addpath(genpath('Auxiliary_Functions'))
addpath(genpath('Estimation_Routines'))
addpath(genpath('Subroutines'))

rng(1);
tic;

% Parallel computing object
num_workers = str2num(getenv('SLURM_CPUS_PER_TASK'));
if ~isempty(num_workers)
    parpool('local', num_workers);
else
    parpool('local');
end
clear num_workers;


%% DECIDE WHICH EXPERIMENT TO RUN

dgp_type = 'G'; % 'MP'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % 'Recursive'; 'IV'; % Either 'ObsShock', 'Recursive', or 'IV'


%% SETTINGS

run(fullfile('Settings', dgp_type));
run(fullfile('Settings', estimand_type));
run(fullfile('Settings', 'shared'));


%% DGP

%----------------------------------------------------------------
% Set up DGP
%----------------------------------------------------------------

% estimate DFM from dataset

DFM_estimate = DFM_est(DF_model.n_fac, DF_model.n_lags_fac, DF_model.n_lags_uar);

% store estimated DFM parameters

DF_model.Phi           = DFM_estimate.Phi;
DF_model.Sigma_eta     = DFM_estimate.Sigma_eta;

DF_model.Lambda        = DFM_estimate.Lambda;
DF_model.delta         = DFM_estimate.delta;
DF_model.sigma_v       = DFM_estimate.sigma_v;

DF_model.variable_name = DFM_estimate.bplabvec_long;

%----------------------------------------------------------------
% Calibrate IV strength and Shock Weight
%----------------------------------------------------------------

% extract factor shock series and external shock series

DFM_estimate.fac_shock                 = DFM_estimate.fac_res / chol(DFM_estimate.Sigma_eta);
external_shock_data                    = readtable(strcat('external_shock_series_', dgp_type, '.csv'));
DFM_estimate.external_shock            = external_shock_data{:,2};
DFM_estimate.external_shock_time_range = [external_shock_data{1,1}, external_shock_data{end,1}, 4];
clear external_shock_data;

% regress external shock series on factor shock series to calibrate

DFM_estimate.calibrate_out       = calibrateIV(DFM_estimate);
DF_model.calibrated_shock_weight = DFM_estimate.calibrate_out.weight;

if DF_model.IV.IV_strength_calibrate==1
    DF_model.IV.alpha = DFM_estimate.calibrate_out.alpha;
    DF_model.IV.sigma_v = DFM_estimate.calibrate_out.sigma_v;
else
    DF_model.IV.alpha = DF_model.IV.manual_alpha;
    DF_model.IV.sigma_v = DF_model.IV.manual_sigma_v;
end

%----------------------------------------------------------------
% Represent as ABCDEF Form
%----------------------------------------------------------------

DF_model.n_s   = size(DF_model.Phi,2);
DF_model.n_eps = size(DF_model.Sigma_eta,2);
DF_model.n_y   = size(DF_model.Lambda,1);
DF_model.n_w   = size(DF_model.delta,1);
DF_model.n_e   = DF_model.n_w * DF_model.n_lags_uar;

DF_model.ABCD  = ABCD_fun_DFM(DF_model);


%% PREPARATION

%----------------------------------------------------------------
% Select Specifications
%----------------------------------------------------------------

settings.specifications = pick_var_fn(DF_model, settings);

%----------------------------------------------------------------
% Results Placeholder
%----------------------------------------------------------------

settings.est.n_methods = length(settings.est.methods_name);
settings.est.full_methods_name = {'svar','svar_corrbias','bvar','lp','lp_penalize','var_avg','svar_iv'};

for i_method = 1:length(settings.est.full_methods_name)
    thisMethod = settings.est.full_methods_name{i_method};
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


%% PRELIMINARY COMPUTATIONS: ESTIMANDS

%----------------------------------------------------------------
% Compute True IRFs in Complete Model
%----------------------------------------------------------------

[DF_model.irf, settings.est.shock_weight] = compute_irfs(DF_model,settings);

%----------------------------------------------------------------
% Compute Degree of Invertibility in Specifications
%----------------------------------------------------------------

DF_model.R0_sq = compute_invert_DFM(DF_model,settings);

%----------------------------------------------------------------
% Compute Persistency of Observables in Specifications
%----------------------------------------------------------------

[DF_model.LRV_Cov_tr_ratio, DF_model.VAR_largest_root, DF_model.frac_coef_for_large_lags] =...
    compute_persist_DFM(DF_model,settings);

%----------------------------------------------------------------
% Compute Target IRF
%----------------------------------------------------------------

switch estimand_type
    
    case 'ObsShock'
        DF_model.target_irf = DF_model.irf(settings.est.IRF_select, ...
            settings.specifications.var_select(:,settings.est.IRF_response_var_pos));
        
    case 'Recursive'
        DF_model.VAR_irf = compute_VARirfs_DFM(DF_model,settings);
        DF_model.target_irf = DF_model.VAR_irf(settings.est.IRF_select, :);
        
    case 'IV'
        DF_model.IV_irf = compute_IVirfs(DF_model,settings);
        DF_model.target_irf = DF_model.IV_irf(settings.est.IRF_select, :);

end

%----------------------------------------------------------------
% Compute Population IV Strengths
%----------------------------------------------------------------

if strcmp(estimand_type, 'IV')
    DF_model.IV_strength = compute_IVstrength_DFM(DF_model, settings);
end


%% MONTE CARLO ANALYSIS

parfor i_MC = 1:settings.simul.n_MC
% for i_MC = 1:settings.simul.n_MC

    if mod(i_MC, 10) == 0
        disp("Monte Carlo:")
        disp(i_MC)
    end

    %----------------------------------------------------------------
    % Generate Data
    %----------------------------------------------------------------

    rng(settings.simul.seed(i_MC));

    data_sim_all = generate_data(DF_model,settings);

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
    
    temp_n_lags_svar = NaN(1,settings.specifications.n_spec);
    temp_n_lags_svar_corrbias = NaN(1,settings.specifications.n_spec);
    temp_n_lags_bvar = NaN(1,settings.specifications.n_spec);
    temp_n_lags_lp = NaN(1,settings.specifications.n_spec);
    temp_n_lags_lp_penalize = NaN(1,settings.specifications.n_spec);
    temp_n_lags_var_avg = NaN(1,settings.specifications.n_spec);
    temp_n_lags_svar_iv = NaN(1,settings.specifications.n_spec);
    
    temp_largest_root_svar = NaN(1,settings.specifications.n_spec);
    temp_LM_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_lambda_lp_penalize = NaN(1,settings.specifications.n_spec);
    temp_weight_var_avg = NaN(2*settings.est.n_lags_max,...
        length(settings.est.average_store_weight),settings.specifications.n_spec);
    temp_F_stat_svar_iv = NaN(1,settings.specifications.n_spec);
    
    %----------------------------------------------------------------
    % Selecting Data
    %----------------------------------------------------------------

    for i_spec = 1:settings.specifications.n_spec
        
        data_sim_select = select_data_fn(data_sim_all,settings,i_spec);
    
        %----------------------------------------------------------------
        % IRF Estimation
        %----------------------------------------------------------------

        % VAR
        
        if any(strcmp(settings.est.methods_name, 'svar'))
            [temp_irf_svar(:,i_spec),temp_n_lags_svar(1,i_spec),temp_largest_root_svar(1,i_spec),temp_LM_stat_svar(1,i_spec)]...
                = SVAR_est(data_sim_select,settings);
        end

        % bias-corrected VAR
        
        if any(strcmp(settings.est.methods_name, 'svar_corrbias'))
            [temp_irf_svar_corrbias(:,i_spec),temp_n_lags_svar_corrbias(1,i_spec)]...
                = SVAR_corr_est(data_sim_select,settings);
        end

        % Bayesian VAR
        
        if any(strcmp(settings.est.methods_name, 'bvar'))
            [temp_irf_bvar(:,i_spec),temp_n_lags_bvar(1,i_spec)]...
                = BVAR_est(data_sim_select,settings);
        end

        % LP

        if any(strcmp(settings.est.methods_name, 'lp'))
            [temp_irf_lp(:,i_spec),temp_n_lags_lp(1,i_spec)]...
                = LP_est(data_sim_select,settings);
        end

        % shrinkage LP

        if any(strcmp(settings.est.methods_name, 'lp_penalize'))
            [temp_irf_lp_penalize(:,i_spec),temp_n_lags_lp_penalize(1,i_spec), temp_lambda_lp_penalize(1,i_spec)]...
                = LP_shrink_est(data_sim_select,settings);
        end

        % VAR model averaging
        
        if any(strcmp(settings.est.methods_name, 'var_avg'))
            [temp_irf_var_avg(:,i_spec),temp_n_lags_var_avg(1,i_spec), temp_weight_var_avg(:,:,i_spec)]...
                = VAR_avg_est(data_sim_select,settings);
        end
        
        % SVAR-IV       

        if any(strcmp(settings.est.methods_name, 'svar_iv'))
            [temp_irf_svar_iv(:,i_spec),temp_n_lags_svar_iv(1,i_spec),temp_F_stat_svar_iv(1,i_spec)]...
                = SVAR_IV_est(data_sim_select,settings);
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
    
    results_n_lags_svar(i_MC,:) = temp_n_lags_svar;
    results_n_lags_svar_corrbias(i_MC,:) = temp_n_lags_svar_corrbias;
    results_n_lags_bvar(i_MC,:) = temp_n_lags_bvar;
    results_n_lags_lp(i_MC,:) = temp_n_lags_lp;
    results_n_lags_lp_penalize(i_MC,:) = temp_n_lags_lp_penalize;
    results_n_lags_var_avg(i_MC,:) = temp_n_lags_var_avg;
    results_n_lags_svar_iv(i_MC,:) = temp_n_lags_svar_iv;
    
    results_largest_root_svar(i_MC,:) = temp_largest_root_svar;
    results_LM_stat_svar(i_MC,:) = temp_LM_stat_svar;
    results_lambda_lp_penalize(i_MC,:) = temp_lambda_lp_penalize;
    results_weight_var_avg(:,:,i_MC,:) = temp_weight_var_avg;
    results_F_stat_svar_iv(i_MC,:) = temp_F_stat_svar_iv;

end

% clear temporary storage

for i_method = 1:length(settings.est.full_methods_name)
    thisMethod = settings.est.full_methods_name{i_method};
    eval(['clear temp_irf_' thisMethod ';']);
    eval(['clear temp_n_lags_' thisMethod ';']);
end
clear temp_largest_root_svar temp_LM_stat_svar temp_lambda_lp_penalize temp_weight_var_avg temp_F_stat_svar_iv
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

if any(strcmp(settings.est.methods_name, 'lp_penalize'))    
    results.lambda.lp_penalize = results_lambda_lp_penalize;
end

if any(strcmp(settings.est.methods_name, 'var_avg'))
    results.weight.var_avg = results_weight_var_avg;
end

if any(strcmp(settings.est.methods_name, 'svar_iv'))
    results.F_stat.svar_iv = results_F_stat_svar_iv;
end

for i_method = 1:length(settings.est.full_methods_name)
    thisMethod = settings.est.full_methods_name{i_method};
    eval(['clear results_irf_' thisMethod ' results_n_lags_' thisMethod ';']);
end

clear results_largest_root_svar results_LM_stat_svar results_lambda_lp_penalize results_weight_var_avg results_F_stat_svar_iv
clear i_method thisMethod

% store IRF only at selected horizons

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
        ' - permute(DF_model.target_irf,[1 3 2])).^2, 2));']);
    
    eval(['results.BIAS2.' thisMethod '= (squeeze(mean(results.irf.' thisMethod ...
        ', 2)) - DF_model.target_irf).^2;']);
    
    eval(['results.VCE.' thisMethod '= squeeze(var(results.irf.' thisMethod ', 0, 2));']);
    
end

clear i_method thisMethod

% export results

mkdir('Results');
save(fullfile('Results', strcat('DFM_', dgp_type, '_', estimand_type)),'DFM_estimate','DF_model','settings','results','-v7.3');
toc;

%% PLOT RESULTS

%----------------------------------------------------------------
% Plot IRFs for Checking
%----------------------------------------------------------------

% for i_method = 1:settings.est.n_methods
%     
%     thisMethod = settings.est.methods_name{i_method};
%     figure(i_method)
%     plot(settings.est.IRF_select, DF_model.target_irf(:,settings.specifications.plot_indx),'Linewidth',5)
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