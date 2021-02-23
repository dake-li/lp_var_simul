%% SPECIFIC SETTINGS FOR DGPs WITH G SHOCKS

% DGP selection

settings.specifications.random_fixed_var      = 12; % always include this variable (= gov't spending) when randomly selecting DGPs
settings.specifications.random_fixed_pos      = 1; % position of fixed variable in each specification

% structural estimands

settings.est.shock_optimize_var_IRF    = 12; % if shock weight is estimated to maximize an IRF, then it is the IRF of this variable in the DFM
settings.est.IRF_response_var_pos      = 2; % interested in IRF of which variable in each specification?
settings.est.IV_est_normalize_var_pos  = 1; % choose IRF normalization variable for all IV methods