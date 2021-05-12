function [IRF,n_lags_est,lambda_opt] = LP_shrink_est(data_sim,settings);
% Function for estimating IRFs using penalized LP

% preparations
lambdaRange = settings.est.lambdaRange;
irfLimitOrder = settings.est.irfLimitOrder;
run('Estimation_Setup'); % common setup for all estimation methods

% data for LP routine
[y, x, w] = LP_gen_data(Y,recursiveShock,responseV,nlags,0); % Warning: here w include both contemperaneous and lagged controls

% settings for LP shrinkage method
H_min = 0; % min horizon in IRF
H_max = IRF_hor - 1; % max horizon in IRF
r = irfLimitOrder + 1; % order of finite difference operator in the penalty term

% leave-one-out cross validation
nT = size(Y,1);
lambdaRange = lambdaRange * nT; % scale the grid of lambda (penalty strength) by nT
lambdaRange = [1e-4, lambdaRange, 1e10]; % allow regular OLS or completely smoothed
rss_cv = locproj_cv(y, x, w, H_min, H_max, r, lambdaRange, settings.est.CV_folds); % cross-validated MSE for each value of lambda
[~,lambda_opt_loc] = min(rss_cv);
lambda_opt = lambdaRange(lambda_opt_loc); % optimally tuned lambda

% re-estimate IRF via penalized LP using the full sample
IRF_resp = locproj(y, x, w, H_min, H_max, r, lambda_opt); % IRF to one unit of shock
IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0); % response of normalization variable estimated by least-squares LP
IRF = IRF_resp / IRF_normalize; % normalize by response of normalization variable

end