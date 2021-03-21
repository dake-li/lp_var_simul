function [IRF,nlags,largest_root,LM_stat,LM_pvalue,Granger_stat,Granger_pvalue] = SVAR_est(data_sim,settings,bias_corrected);
% Function for estimating IRFs using least-squares VAR or bias-corrected VAR

% preparations

res_autocorr_nlags = settings.est.res_autocorr_nlags;
run('Estimation_Setup'); % common setup for all estimation methods

% estimate VAR

if bias_corrected == 0
    [Bc,By,Sigma,Sxx,Res,Beta] = VAR(Y,nlags); % no bias correction
else
    [Bc,By,Sigma] = VAR_CorrectBias(Y,nlags); % with bias correction
end

G = chol(Sigma, 'lower'); % Warning: correspond to matrix C in our paper
ShockVector = G(:,recursiveShock);

% estimate IRF

IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1); % IRF to one unit of shock
IRF = IRF(responseV,:) / IRF(normalizeV,1); % normalize by response of normalization variable
IRF = IRF';

% when largest_root, LM_stat, Granger_stat are computed

if nargout > 2
    
    if bias_corrected == 1 % not compute diagonosis for bias-corrected VAR
        
        largest_root = NaN;
        LM_stat = NaN;
        LM_pvalue = NaN;
        Granger_stat = NaN;
        Granger_pvalue = NaN;
        
    else % compute diagonosis for simple VAR
        
        % estimate largest root in VAR

        nv = size(Y,2);
        companion_form = diag(ones(1, nv*(nlags-1)), -nv);
        companion_form(1:nv,:) = reshape(By,[nv, nv*nlags]);
        largest_root = max(abs(eig(companion_form)));

        % estimate LM-stat to examine VAR(p) fit (p is fixed or chosen by AIC)

        nT = size(Y,1);
        Y_lag = lagmatrix(Y,1:nlags); % lagged Y as explanatory variables
        Y_lag = Y_lag((nlags+res_autocorr_nlags+1):end,:);
        Res_lag = lagmatrix(Res,1:res_autocorr_nlags); % lagged residual
        Res_lag = Res_lag((res_autocorr_nlags+1):end, :);
        X_auxiliary = [ones(nT-nlags-res_autocorr_nlags,1), Y_lag, Res_lag]; % auxiliary controls
        Res_current = Res((res_autocorr_nlags+1):end, :); % current residual
        [~,~,~,Res_aux] = LS(Res_current, X_auxiliary); % auxiliary regression
        Sigma_original = cov(Res_current,1); % res var-cov matrix in the original regression
        Sigma_auxiliary = cov(Res_aux,1); % res var-cov matrix in the auxiliary regression
        LM_stat = (nT-nlags-res_autocorr_nlags - nv*nlags - 1 - nv*res_autocorr_nlags - 0.5) *...
            log(det(Sigma_original) / det(Sigma_auxiliary));
        LM_pvalue = chi2cdf(LM_stat, nv^2 * res_autocorr_nlags, 'upper');

        % test if IV granger-causes endogenous variables
        
        if with_IV == 0
            
            % return NaN if there is no IV
            
            Granger_stat = NaN;
            Granger_pvalue = NaN;
            
        else
            
            % return Wald stat for all the VAR coefs of IV being jointly
            % zero in all the equations of w_t
            
            nv = size(Y,2);
            nT = size(Y,1);
            beta = Beta'; % dim: nv * (nv * nlags + 1)
            beta_IV_loc = zeros(size(beta));
            beta_IV_loc(2:nv, 1 + (0:(nlags-1))*nv + 1) = 1; % indicator of the IV coefs
            beta_IV_loc = logical(beta_IV_loc);
            beta_IV = beta(beta_IV_loc); % extract IV coefs
            Omega = kron(inv(Sxx), Sigma);
            Omega_IV = Omega(beta_IV_loc(:), beta_IV_loc(:)); % asymptotic variance matrix of IV coefs
            Granger_stat = (nT-nlags)* beta_IV' / Omega_IV * beta_IV; % Wald stat to jointly test IV coefs
            Granger_pvalue = chi2cdf(Granger_stat, (nv-1) * nlags, 'upper');
        end
    
    end
    
end

end