function [IRF,nlags] = SVAR_corr_est(data_sim,settings);

% preparations

run('Estimation_Setup');

% estimate VAR

[~,ByCorrect,Sigma,~] = VAR_CorrectBias(Y,nlags);
G = chol(Sigma, 'lower');
ShockVector = G(:,recursiveShock);
IRF = IRF_SVAR(ByCorrect,ShockVector,IRF_hor - 1);
IRF = IRF(responseV,:) / IRF(normalizeV,1);
IRF = IRF';

end