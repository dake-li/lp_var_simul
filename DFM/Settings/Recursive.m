%% SPECIFIC SETTINGS FOR RECURSIVE SHOCK IDENTIFICATION

% indicate position of recursive shock of interest

settings.est.recursive_shock_pos = settings.specifications.random_fixed_pos; % which is recursively defined shock?

% indicate estimand

settings.est.with_shock      = 0; % shock is observed and ordered first in data?
settings.est.recursive_shock = 1; % use recursive shock
settings.est.with_IV         = 0; % IV is ordered first in data or used in IV method?