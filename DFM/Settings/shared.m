%% ENCOMPASSING DFM

% take DFM dimensions from Stock-Watson (2016)

DF_model.n_y        = 207; % number of observables
DF_model.reorder    = [1:76, 87:94, 77:86, 95:171, 181:195, 172:180, 196:207]; % index to reorder data to match variable list in Stock-Watson (2016)

DF_model.levels     = 1; % =1: variables in levels, 0=: differenced
DF_model.n_fac      = 6; % number of factors
DF_model.n_lags_fac = 2; % lag order of factors
DF_model.n_lags_uar = 2; % lag order of measurement error

% scale up calibrated persistency of factors

DF_model.fac_persist.scale = 0; % scale calibrated persistence of factors?
DF_model.fac_persist.target_largest_root = 0.95; % match the targeted largest root in factor equation
DF_model.fac_persist.keep_fac_variation = 1; % adjust factor innovation scale to keep factor variation?

%% PREPARATIONS FOR STRUCTURAL ESTIMANDS

% selection of DGPs from encompassing model

settings.specifications.random_select         = 1; % randomly select variables from DFM list?
settings.specifications.random_n_spec         = 7; % number of random specifications
settings.specifications.random_n_var          = 5; % number of variables in each random specification
settings.specifications.random_category_range = [1 20; 21 31; 32 76; 77 86; 87 94; 95 131; 132 141;...
                                                 142 159; 160 171; 172 180; 181 207]; % ranges for Stock-Watson variable categories (see their Table 1)
settings.specifications.random_category_setup = {[1,2,3], 6}; % at least draw one from certain categories
settings.specifications.random_from_key_series = 0; % randomly select from some key series in DFM list?
settings.specifications.random_exhaust_key_series = 0; % exhaust all the combinations from key series, or have only random_n_spec draws?
settings.specifications.random_key_series = [1,2,6,12,56,95,97,121,132,142,147,151,172,181,193,196,202]; % list of key series

settings.specifications.manual_var_select     = [1 142; 1 97]; % if manual selection, then these sets of variables will be selected

% preliminary structural shock settings for observed shock and IV

settings.est.estimate_shock_weight    = 1; % do we estimate the loading of the structural shock on the reduced-form shocks?
settings.est.shock_weight_calibrate   = 0; % if shock weight is estimated: 0 to maximize spending/rate IRF, 1 to calibrate to external shocks

settings.est.manual_shock_pos         = 1; % if loading is not estimated, we just mechanically choose the xth shock in the factor VAR

% IRFs of interest

settings.est.IRF_hor              = 21; % maximal horizon (include contemporary)
settings.est.IRF_select           = 1:settings.est.IRF_hor; % which IRFs to study

% compute R0_sq using VMA representation

settings.est.VMA_nlags = 50; % number of lags in VMA representation

% compute largest root and VAR(p) fit in population using truncated infinite-order VAR 

settings.est.VAR_infinity_truncate = 50; % truncation horizon for VAR(infty)

% number of Monte Carlo draws

settings.simul.n_MC    = 4; % number of Monte Carlo reps
settings.simul.seed    = (1:settings.simul.n_MC)*10 + randi([0,9],1,settings.simul.n_MC); % random seed for each Monte Carlo

% sample settings

settings.simul.T      = 200; % time periods for each simulation
settings.simul.T_burn = 100; % burn-in

%% ESTIMATION SETTINGS

% set of estimation methods

settings.est.methods_name    = {'svar','svar_corrbias','bvar','lp','lp_penalize','var_avg'};

% lag specification

settings.est.est_n_lag      = isnan(lag_type); % do we use the estimated lag order, or just set it?
settings.est.est_n_lag_BIC  = 0; % lag order estimation: use BIC? default is AIC
settings.est.n_lags_fix     = lag_type; % default number of lags if not estimated
settings.est.n_lags_max     = 20; % maximal lag length for info criteria

settings.est.res_autocorr_nlags = 2; % check autocorr of VAR(p) residuals up to this order

% BVAR prior

settings.est.prior.towards_random_walk = DF_model.levels; % prior shrinking towards random walk? otherwise towards zero
settings.est.prior.tight_overall       = 0.04;
settings.est.prior.tight_nonown_lag    = 0.25;
settings.est.prior.decay_power         = 2;
settings.est.prior.tight_exogenous     = 1e5;

% BVAR posterior draws

settings.est.posterior_ndraw = 100; % number of posterior draws (if set to 0, only use posterior mean of VAR coefficient to compute posterior mean of IRF)

% LP smoothing

settings.est.lambdaRange   = [0.001:0.005:0.021, 0.05:0.1:1.05, 2:1:19, 20:20:100, 200:200:2000]; % cross validation grid, scaled up by T
settings.est.irfLimitOrder = 2; % shrink towards polynomial of that order
settings.est.CV_folds      = 5; % Number of folds used for cross validation

% VAR model averaging

settings.est.average_store_weight       = [2, 5, 9, 15, 21]; % store model weights at which horizon
settings.est.average_store_submodel_irf = 0; % store IRF of each submodel? Only store if want to comput oracle weight
settings.est.average_max_lags           = 1; % include lags up to n_lags_max? otherwise up to estimated lags
settings.est.average_options            = optimoptions('quadprog','Display','off');
