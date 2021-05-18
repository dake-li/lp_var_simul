%% SPECIFIC SETTINGS FOR IV IDENTIFICATION

% IV equation

DF_model.IV.manual_rho     = 0.1; % manually set up baseline IV persistence
DF_model.IV.manual_alpha   = 1; % manually set up IV shock coefficient
DF_model.IV.manual_sigma_v = 1; % manually set up IV noise

settings.est.IV.IV_persistence_calibrate = 1; % use calibrated baseline IV persistence
settings.est.IV.IV_persistence_scale = [0.5 1 2]; % scale up or down baseline persistence (use 1 by default)
settings.est.IV.IV_strength_calibrate = 1; % use calibrated IV strength (alpha = 1, sigma_v will be changed)

% estimation methods

settings.est.methods_name    = [settings.est.methods_name {'svar_iv'}]; % add SVAR-IV to the estimators

% indicate estimand

settings.est.with_shock      = 0; % shock is observed and ordered first in data?
settings.est.recursive_shock = 0; % use recursive shock
settings.est.with_IV         = 1; % IV is ordered first in data or used in IV method?