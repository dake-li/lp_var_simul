%% SPECIFIC SETTINGS FOR OBSERVED SHOCK IDENTIFICATION

% indicate estimand

settings.est.with_shock      = 1; % shock is observed and ordered first in data?
settings.est.recursive_shock = 0; % use recursive shock
settings.est.with_IV         = 0; % IV is ordered first in data or used in IV method?

% indicate normalization scheme
settings.est.normalize_with_shock_std_dev = 1; % normalize IRF with one unit of shock (shock std-dev is 1)? Otherwise using normalization variable