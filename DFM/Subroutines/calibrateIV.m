function out = calibrateIV(DFM_estimate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% load factor shock
factor_shock = DFM_estimate.fac_shock;
factor_shock_time_range = DFM_estimate.factor_shock_time_range;

% cut NaN at the beginning of data
number_of_NaN = sum(isnan(factor_shock(:,1)));
factor_shock = factor_shock((number_of_NaN+1):end, :);
factor_shock_time_range(1) = factor_shock_time_range(1) + number_of_NaN / factor_shock_time_range(3);

% load external shock
external_shock = DFM_estimate.external_shock;
external_shock_time_range = DFM_estimate.external_shock_time_range;

% pick overlapping time range
overlap_time_range = [NaN, NaN, factor_shock_time_range(3)];
overlap_time_range(1) = max(factor_shock_time_range(1), external_shock_time_range(1));
overlap_time_range(2) = min(factor_shock_time_range(2), external_shock_time_range(2));

% pick rows of data
factor_shock_index = [NaN, NaN];
factor_shock_index(1) = (overlap_time_range(1) - factor_shock_time_range(1)) * factor_shock_time_range(3) + 1;
factor_shock_index(2) = (overlap_time_range(2) - factor_shock_time_range(1)) * factor_shock_time_range(3) + 1;
external_shock_index = [NaN, NaN];
external_shock_index(1) = (overlap_time_range(1) - external_shock_time_range(1)) * external_shock_time_range(3) + 1;
external_shock_index(2) = (overlap_time_range(2) - external_shock_time_range(1)) * external_shock_time_range(3) + 1;

% construct regressors
X = factor_shock(factor_shock_index(1):factor_shock_index(2),:);
X = [ones(size(X,1),1),X];
Y = external_shock(external_shock_index(1):external_shock_index(2),1);

% regress external shock on factor shocks
[Beta,Sigma,~,~] = LS(Y,X);
R2 = 1 - Sigma / var(Y); % adjusted R2
weight = Beta(2:end);
alpha = sqrt(sum(weight.^2)); % IV shock coefficient
weight = weight / sqrt(sum(weight.^2)); % weight on factor shock
sigma_v = sqrt(Sigma); % IV noise
nobs = size(Y,1); % number of observation

% pack up result
out.weight = weight;
out.alpha = alpha;
out.sigma_v = sigma_v;
out.R2 = R2;
out.nobs = nobs;

end

