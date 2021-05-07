function [stat, pvalue, irf_LP, varcov_VAR, varcov_LP] = IRF_Hausman(Y,respV,irf_VAR,By,Sigma,Sxx)

    % Hausman test of null hypothesis: SVAR IRF = LP IRF
    
    % Test only includes horizons greater than "nlags"
    % Var-cov matrix of LP IRF computed under the null of well-specified VAR(p) model
    % Code only applies to the case with an observed shock (first column of Y)
    
    
    % Dimensions
    nv = size(By,1);
    nlags = size(By,3);
    nhorizons = length(irf_VAR)-1;
    
    % Compute LP IRF
    irf_LP = IRF_LP(Y,1,respV,nlags,nhorizons)';
    
    % Compute various IRFs based on VAR
    G = chol(Sigma, 'lower'); % Warning: corresponds to matrix B in our paper
    gamma = G(:,1)/G(1,1);
    irf_gamma_VAR = IRF_SVAR(By,gamma,nhorizons); % Normalized IRFs of all variables wrt. first orthogonalized shock
    irf_respV_redf_VAR = nan(nv,1+nhorizons); % Will contain reduced-form IRFs of response variable
    the_eye = eye(nv);
    for i=1:nv
        the_irf_redf = IRF_SVAR(By,the_eye(:,i),nhorizons);
        irf_respV_redf_VAR(i,:) = the_irf_redf(respV,:);
    end
    irf_respV_VAR = G'*irf_respV_redf_VAR; % IRFs of response variable wrt. all orthogonalized shocks
    
    % Build LP residual MA structure and VAR IRF Jacobian
    jacob_VAR = zeros(nv,nv*nv*nlags,1+nhorizons); % Will contain Jacobian of VAR IRF wrt. By (VAR coefficients)
    ma_coef_LP = zeros(nv*(1+nhorizons),1+nhorizons); % Will contain MA loadings of LP residuals on shocks epsilon_t, epsilon_{t+1}, ...
    
    for h=0:nhorizons
        
        % MA structure for LP residual with response variable y_{t+h}
        ma_coef_LP(1:nv*(1+h),1+h) = reshape(irf_respV_VAR(:,1+h:-1:1),1,[]); % LP residual MA coef's = orthogonalized IRFs
        ma_coef_LP(1,1+h) = 0; % epsilon_{1,t} does not appear in residual since we regress on it
        
        % Compute VAR Jacobian via chain rule applied to IRF recursion
        aux = zeros(nv,nv*nv*nlags);
        for j=1:min(h,nlags)
            aux = aux + By(:,:,j)*jacob_VAR(:,:,1+h-j);
            aux(:,nv*nv*(j-1)+1:nv*nv*j) = aux(:,nv*nv*(j-1)+1:nv*nv*j) ...
                                           + kron(irf_gamma_VAR(:,1+h-j)',eye(nv));
        end
        jacob_VAR(:,:,1+h) = aux;
        
    end
    
    % LP var-cov matrix implied by estimated VAR(p) model
    varcov_LP = ma_coef_LP'*ma_coef_LP/(Y(nlags+nhorizons+1:end,1)'*Y(nlags+nhorizons+1:end,1));
    
    % VAR var-cov matrix
    Omega = kron(inv(Sxx), Sigma)/(size(Y,1)-nlags); % Var-cov matrix of vec([Bc By])
    varcov_LP_impact_other = (G(2:end,2:end)*G(2:end,2:end)')/(Y(nlags+1:end,1)'*Y(nlags+1:end,1));
        % Var-cov matrix of LP impact impulse responses for variables 2, 3, ...
    jacob_respV_VAR = permute(jacob_VAR(respV,:,:), [2 3 1]); % Jacobian for response variable
    varcov_VAR = jacob_respV_VAR'*Omega(nv+1:end,nv+1:end)*jacob_respV_VAR ...
                 + irf_respV_redf_VAR(2:end,:)'*varcov_LP_impact_other*irf_respV_redf_VAR(2:end,:);
    
    % Hausman test statistic and p-value
    irf_diff = irf_VAR-irf_LP; % IRF difference
    varcov_diff = varcov_LP-varcov_VAR; % Var-cov matrix of IRF difference (under the null)
    stat = irf_diff(nlags+2:end)'*(varcov_diff(nlags+2:end,nlags+2:end)\irf_diff(nlags+2:end));
    df = nhorizons-nlags;
    pvalue = 1-chi2cdf(stat,df);
    

end