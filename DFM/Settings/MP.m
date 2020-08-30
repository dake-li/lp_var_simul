% Settings: dgp_type = MP


%% Specification

settings.specifications.random_fixed_var      = 142; % always include which variable when random select
settings.specifications.random_fixed_pos      = 5; % position of fixed variable in each specification


%% Estimation

settings.est.IRF_response_var_pos = 1; % interested in IRF of which variable in each specification?

settings.est.IV_est_normalize_var_pos = 5; % choose IRF normalization variable in all IV methods