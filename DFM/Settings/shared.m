%% DGP

% Stock-Watson DFM dimensions

DF_model.n_y        = 207; % number of observables

DF_model.n_fac      = 6; % number of factor
DF_model.n_lags_fac = 2; % lag order of factor
DF_model.n_lags_uar = 2; % lag order of measurement error


%% Experiment Specification

% variable selection

settings.specifications.manual_var_select     = [1 142; 1 97]; % manually select specifications
settings.specifications.random_select         = 1; % randomly select?
settings.specifications.random_n_spec         = 2; % number of random specifications
settings.specifications.random_n_var          = 5; % number of variables in each random specification
settings.specifications.random_category_range = [1 20; 21 31; 32 76; 77 86; 87 94; 95 131; 132 141;...
                                                 142 159; 160 171; 172 180; 181 207];
                                                   % category ranges
settings.specifications.random_category_setup = {[1,2,3], 6}; % at least draw one from certain categories
settings.specifications.plot_indx             = 1; % plot the only specification

% shock position

settings.est.manual_shock_pos         = 1; % manually choose which shock to be our true structural shock?
settings.est.estimate_shock_weight    = 1; % automatically estimate shock weights for true shock? 
settings.est.shock_weight_calibrate = 0; % when estimate shock weight, use calibrated result? or optimize targeted IRF

% IRFs of interest

settings.est.IRF_hor              = 20; % maximal horizon (include contemporary)
settings.est.IRF_select           = 1:20; % which IRFs to summarize

% compute R0_sq using VMA representation

settings.est.VMA_nlags = 50;

% compute largest root and VAR(p) fit in population using truncated infinite-order VAR 

settings.est.VAR_infinity_truncate = 50; 

% number of Monte Carlo draws

settings.simul.n_MC    = 3; % number of Monte Carlo reps
settings.simul.seed    = (1:settings.simul.n_MC)*10 + randi([0,9],1,settings.simul.n_MC); % random seed for each Monte Carlo

% simulation details

settings.simul.T      = 200; % time periods for each simulation
settings.simul.T_burn = 100; % burn-in


%% Estimation Settings

% choose estimand

settings.est.methods_name    = {'svar','svar_corrbias','bvar','lp','lp_penalize','var_avg'}; % choose estimands (may be expanded later)

% lag specification

settings.est.est_n_lag      = isnan(lag_type); % estimate number of lags?
settings.est.est_n_lag_BIC  = 0; % use BIC? otherwise use AIC
settings.est.n_lags_fix     = lag_type; % default number of lags if not estimated
settings.est.n_lags_max     = 20; % maximal lag length for info criteria
settings.est.res_autocorr_nlags = 1; % check autocorr of VAR(p) residuals up to this order

% BVAR prior

settings.est.prior.tight_overall      = 0.04;
settings.est.prior.tight_nonown_lag   = 0.25;
settings.est.prior.decay_power        = 2;
settings.est.prior.tight_exogenous    = 1e5;

% LP smoothing

settings.est.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, 2:1:19, 20:20:100, 200:200:2000]; % cross validation grid, scaled up by T
settings.est.irfLimitOrder = 2; % shrink towards polynomial of that order
settings.est.CV_folds      = 5; % Number of folds used for cross validation

% VAR model averaging

settings.est.average_store_weight = [2, 11, 20]; % store model weights at which horizon
settings.est.average_max_lags     = 1; % include lags up to n_lags_max? otherwise up to estimated lags
settings.est.average_options      = optimoptions('quadprog','Display','off');
