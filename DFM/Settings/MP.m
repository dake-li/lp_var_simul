%% SPECIFIC SETTINGS FOR DGPs WITH MP SHOCKS

% DGP selection

settings.specifications.random_fixed_var      = 142; % always include this variable (= fed funds rate) when randomly selecting DGPs
settings.specifications.random_fixed_pos      = settings.specifications.random_n_var; % position of fixed variable in each specification (= end)

% structural estimands

settings.est.shock_optimize_var_IRF    = 142; % if shock weight is estimated to maximize an IRF, then it is the IRF of this variable in the DFM
settings.est.IRF_response_var_pos      = 1; % interested in IRF of which variable in each specification?
settings.est.IV_est_normalize_var_pos  = 5; % choose IRF normalization variable for all IV methods