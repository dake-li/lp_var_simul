function [IRF,nlags] = BVAR_est(data_sim,settings);

% preparations

prior = settings.est.prior;
run('Estimation_Setup');

% estimate VAR

[~,By,Sigma,~] = BVAR(Y,nlags,prior);
G = chol(Sigma, 'lower');
ShockVector = G(:,recursiveShock);
IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1);
IRF = IRF(responseV,:) / IRF(normalizeV,1);
IRF = IRF';

end