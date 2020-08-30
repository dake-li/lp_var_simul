% Settings: dgp_type = G


%% Specification

settings.specifications.random_fixed_var      = 12; % always include which variable when random select
settings.specifications.random_fixed_pos      = 1; % position of fixed variable in each specification


%% Estimation

settings.est.IRF_response_var_pos = 2; % interested in IRF of which variable in each specification?

settings.est.IV_est_normalize_var_pos = 1; % choose IRF normalization variable in all IV methods