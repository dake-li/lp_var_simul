function [IRF,nlags,lambda_opt] = LP_shrink_est(data_sim,settings);

% preparations

lambda = settings.est.lambda;
lambdaRange = settings.est.lambdaRange;
irfLimitOrder = settings.est.irfLimitOrder;
run('Estimation_Setup');

% data for LP routine
[y, x, w] = LP_gen_data(Y,recursiveShock,responseV,nlags,0);

% settings for LP shrinkage method
H_min = 0;
H_max = IRF_hor - 1;
r = irfLimitOrder + 1;

% leave-one-out cross validation
nT = size(Y,1);
lambdaRange = lambdaRange * nT;
lambdaRange = [1e-4, lambdaRange, 1e10]; % allow regular OLS or completely smoothed
rss_cv = locproj_cv(y, x, w, H_min, H_max, r, lambdaRange, settings.est.CV_folds);
[~,lambda_opt_loc] = min(rss_cv);
lambda_opt = lambdaRange(lambda_opt_loc);

% re-estimate
IRF_resp = locproj(y, x, w, H_min, H_max, r, lambda_opt);
IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0);
IRF = IRF_resp / IRF_normalize;

end