function [IRF,nlags] = BVAR_est(data_sim,settings);
% Function for estimating IRFs using a Bayesian VAR approach

% preparations

prior = settings.est.prior;
run('Estimation_Setup'); % common setup for all estimation methods

% estimate Bayesian VAR

[~,By,Sigma,~] = BVAR(Y,nlags,prior);
G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
ShockVector = G(:,recursiveShock);

% estimate IRF

IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1); % IRF to one unit of shock
IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable
IRF = IRF';

end