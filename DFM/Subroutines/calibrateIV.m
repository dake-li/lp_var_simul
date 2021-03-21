function out = calibrateIV(DFM_estimate)
% Function for calibrating IV strength in the IV DGP
    % IV DGP: z_t = \rho z_{t-1} + \alpha (shock_weight' * \epsilon_t) + \nu_t
    % where \alpha = 1, Var(\nu_t) is tuned to match the signal-noise ratio in
    % the calibration regression
    
    % Calibration regression : z_t = c + \sum_{l=-p}^p (a_l' * \epsilon_{t-l}) + residuals_t
    % where z_t is the observed external IV, \epsilon is the fitted residual in DFM
    %       R^2 represents the degree of recoverability of the IV shock using the factors
    
    % Match Var(\nu_t) = 1/R^2 - 1

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

% iterate thru different number of leads and lags
max_num = 10; % max number of leads and lags
BIC_list = NaN(1, max_num+1);
Sigma_list = NaN(1, max_num+1);

for i_num = 0:max_num

    % construct regressors
    Y = external_shock(external_shock_index(1):external_shock_index(2),1);
    nobs = size(Y,1); % number of observations
    X = factor_shock(factor_shock_index(1):factor_shock_index(2),:);
    X = lagmatrix(X, (-i_num):i_num);
    X = [ones(size(X,1),1),X];
    Y = Y((max_num+1):(nobs-max_num),1);
    X = X((max_num+1):(nobs-max_num),:);


    % regress external shock on factor shocks
    [~,Sigma,~,~] = LS(Y,X);
    
    % BIC
    BIC = (nobs- 2 * max_num) * log(Sigma) + (2 * max_num + 1 + 1) * log(nobs- 2 * max_num);
    BIC_list(i_num+1) = BIC;
    Sigma_list(i_num+1) = Sigma;

end

% BIC chooses the optimal number of leads and lags
[~,optimal_num] = min(BIC_list);
optimal_Sigma = Sigma_list(optimal_num);
optimal_num = optimal_num - 1;
R2 = 1 - optimal_Sigma / var(Y); % adjusted R2
alpha = 1; % IV shock coefficient
sigma_v = sqrt((1-R2)/R2) * alpha; % IV noise

% pack up result
out.rho = []; % placeholders for calibrated rho
out.weight = []; % placeholders for calibrated weight
out.alpha = alpha;
out.sigma_v = sigma_v;
out.R2 = R2;
out.nobs = nobs;
out.optimal_num_lead_lag = optimal_num;

end

