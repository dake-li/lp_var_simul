%% LP vs VAR: DFM SIMULATION STUDY
% Dake Li and Christian Wolf

%% HOUSEKEEPING

clc
clear all
close all

addpath(genpath(fullfile('..', '../Auxiliary_Functions')))
addpath(genpath(fullfile('..', '../Estimation_Routines')))
addpath(genpath('../Subroutines'))

rng(1, 'twister');

%% DECIDE WHICH EXPERIMENT TO RUN

% manually set up experiment

spec_id = 1; % specification choice set id (i.e. the seed for random draws of specifications)
dgp_type = 'G'; % 'MP'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % 'Recursive'; 'IV'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)

%% SETTINGS

% Apply shared settings as well as settings specific to DGP and estimand type

run(fullfile('../Settings', 'shared'));
run(fullfile('../Settings', dgp_type));
run(fullfile('../Settings', estimand_type));

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

if strcmp(estimand_type, 'IV')
    if DF_model.IV.IV_strength_calibrate==1
        DF_model.IV.alpha = DFM_estimate.calibrate_out.alpha;
        DF_model.IV.sigma_v = DFM_estimate.calibrate_out.sigma_v;
    else
        DF_model.IV.alpha = DF_model.IV.manual_alpha;
        DF_model.IV.sigma_v = DF_model.IV.manual_sigma_v;
    end
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

settings.specifications = pick_var_fn(DF_model, settings, spec_id);

%----------------------------------------------------------------
% Results Placeholder
%----------------------------------------------------------------

settings.est.n_methods = length(settings.est.methods_name);

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

%% PLOT RESULTS

if strcmp(dgp_type,'G')

spec_select = [1 2 3];

figure
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
for i_spec_indx = 1:length(spec_select)
    i_spec = spec_select(i_spec_indx);
    plot(settings.est.IRF_select, DF_model.target_irf(:,i_spec) ./ max(abs(DF_model.target_irf(:,i_spec))),'Linewidth',3.5)
    hold on
end
hold off
set(gcf,'color','w')
title('Observed G Shock','interpreter','latex','fontsize',20)
xlabel('Horizon','interpreter','latex','FontSize',20)
% ylabel('\% deviation','interpreter','latex','FontSize',20)
xlim([1 20])
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('DFM_illustration_G','-depsc');

elseif strcmp(dgp_type,'MP')
    
spec_select = [1 2 3];

figure
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
for i_spec_indx = 1:length(spec_select)
    i_spec = spec_select(i_spec_indx);
    plot(settings.est.IRF_select, DF_model.target_irf(:,i_spec) ./ max(abs(DF_model.target_irf(:,i_spec))),'Linewidth',3.5)
    hold on
end
hold off
set(gcf,'color','w')
title('Recursive MP Shock','interpreter','latex','fontsize',20)
xlabel('Horizon','interpreter','latex','FontSize',20)
% ylabel('\% deviation','interpreter','latex','FontSize',20)
xlim([1 20])
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('DFM_illustration_MP','-depsc');

end