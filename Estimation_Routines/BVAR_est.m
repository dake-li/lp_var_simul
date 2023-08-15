function [IRF,n_lags_est,GLP_hyper_est] = BVAR_est(data_sim,settings);
% Function for estimating IRFs using a Bayesian VAR approach

% preparations

prior = settings.est.prior;
ndraw = settings.est.posterior_ndraw;
run('Estimation_Setup'); % common setup for all estimation methods

% estimate Bayesian VAR and get posterior mean and variance for VAR coef.

if settings.est.bvar_glp == 1 % If Giannone, Lenza & Primiceri procedure

    n_Y = size(Y,2);

    % Estimate BVAR
    if settings.est.prior.towards_random_walk % Random walk prior
        r = bvarGLP(Y,nlags,'MNpsi',0,'Fcast',0);
    else % White noise prior
        r = bvarGLP(Y,nlags,'MNpsi',0,'Fcast',0,'sur',0,'noc',0,'posi',1:n_Y);
    end

    % Posterior draws of parameters
    if ndraw == 0 % only use posterior means
        beta_draws = r.postmax.betahat;
        sigma_draws = r.postmax.sigmahat;
    else % draw from posterior
        beta_draws = nan(1+n_Y*nlags,n_Y,ndraw);
        sigma_draws = nan(n_Y,n_Y,ndraw);
        for j=1:ndraw
            [beta_draws(:,:,j),sigma_draws(:,:,j)] = post_draw(r.postmax.betahat,r.postmax.Sinv,r.postmax.cholZZinv,r.postmax.T);
        end
    end

    % IRFs
    IRF_draws = nan(n_Y,IRF_hor,size(beta_draws,3));
    for j=1:size(IRF_draws,3)
        G = chol(sigma_draws(:,:,j),'lower');
        ShockVector = G(:,recursiveShock);
        IRF_draws(:,:,j) = IRF_SVAR(reshape(beta_draws(2:end,:,j)',n_Y,n_Y,nlags),ShockVector,IRF_hor - 1);
    end

    % hyperparameters

    GLP_hyper_est = [r.postmax.lambda;r.postmax.theta;r.postmax.miu];

else % Otherwise, basic MN prior with fixed hyper-parameters

    if ndraw == 0 % only use posterior mean of VAR coef
        [~,By,Sigma] = BVAR(Y,nlags,prior);
    else % use posterior mean and variance of VAR coef
        [~,~,Sigma,~,VAR_coef_post_mean,VAR_coef_post_vce_inv] = BVAR(Y,nlags,prior);
    end
    G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
    ShockVector = G(:,recursiveShock);
    
    % estimate IRF
    
    if ndraw == 0 % only use posterior mean of VAR coef to compute IRF
        IRF_draws = IRF_SVAR(By,ShockVector,IRF_hor - 1); % IRF to one unit of shock
    else % compute posterior mean of IRF based on posterior draws
        IRF_draws = IRF_BVAR(VAR_coef_post_mean,VAR_coef_post_vce_inv,ShockVector,IRF_hor - 1,ndraw); % IRF draws
    end

    GLP_hyper_est = NaN;

end

IRF = mean(IRF_draws,3); % posterior mean of IRF draws (or just single IRF if ndraw=0)

IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable
IRF = IRF';

end