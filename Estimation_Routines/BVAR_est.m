function [IRF,n_lags_est] = BVAR_est(data_sim,settings);
% Function for estimating IRFs using a Bayesian VAR approach

% preparations

prior = settings.est.prior;
ndraw = settings.est.posterior_ndraw;
run('Estimation_Setup'); % common setup for all estimation methods

% estimate Bayesian VAR and get posterior mean and variance for VAR coef.

if ndraw == 0 % only use posterior mean of VAR coef
    [~,By,Sigma] = BVAR(Y,nlags,prior);
else % use posterior mean and variance of VAR coef
    [~,~,Sigma,~,VAR_coef_post_mean,VAR_coef_post_vce_inv] = BVAR(Y,nlags,prior);
end
G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
ShockVector = G(:,recursiveShock);

% estimate IRF

if ndraw == 0 % only use posterior mean of VAR coef to compute IRF
    IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1); % IRF to one unit of shock
else % compute posterior mean of IRF based on posterior draws
    IRF_draws = IRF_BVAR(VAR_coef_post_mean,VAR_coef_post_vce_inv,ShockVector,IRF_hor - 1,ndraw); % IRF draws
    IRF = mean(IRF_draws,3); % posterior mean of IRF draws
end

IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable
IRF = IRF';

end