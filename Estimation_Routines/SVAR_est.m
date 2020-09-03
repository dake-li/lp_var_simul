function [IRF,nlags,largest_root,LM_stat] = SVAR_est(data_sim,settings,bias_corrected);

% preparations

res_autocorr_nlags = settings.est.res_autocorr_nlags;
run('Estimation_Setup');

% estimate VAR

if bias_corrected == 0
    [~,By,Sigma,~,Res] = VAR(Y,nlags); % no bias correction
else
    [~,By,Sigma,~,Res] = VAR_CorrectBias(Y,nlags); % with bias correction
end

G = chol(Sigma, 'lower');
ShockVector = G(:,recursiveShock);

% estimate IRF

IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1);
IRF = IRF(responseV,:) / IRF(normalizeV,1);
IRF = IRF';

% when largest_root and LM_stat are computed

if nargout > 2
    
    % estimate largest root in VAR

    nv = size(Y,2);
    companion_form = diag(ones(1, nv*(nlags-1)), -nv);
    companion_form(1:nv,:) = reshape(By,[nv, nv*nlags]);
    largest_root = max(abs(eig(companion_form)));

    % estimate LM-stat to examine VAR(p) fit

    nT = size(Y,1);
    Y_lag = lagmatrix(Y,1:nlags); % lagged Y as explanatory variables
    Y_lag = Y_lag((nlags+res_autocorr_nlags+1):end,:);
    Res_lag = lagmatrix(Res,1:res_autocorr_nlags); % lagged residual
    Res_lag = Res_lag((res_autocorr_nlags+1):end, :);
    X_auxiliary = [ones(nT-nlags-res_autocorr_nlags,1), Y_lag, Res_lag];
    Res_current = Res((res_autocorr_nlags+1):end, :); % current residual
    [~,~,~,Res_aux] = LS(Res_current, X_auxiliary);
    Sigma_original = cov(Res_current,1);
    Sigma_auxiliary = cov(Res_aux,1);
    LM_stat = (nT-nlags-res_autocorr_nlags - nv*nlags - 1 - nv*res_autocorr_nlags - 0.5) *...
        log(det(Sigma_original) / det(Sigma_auxiliary));
    
end

end