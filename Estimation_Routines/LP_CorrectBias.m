function IRF_corr = LP_CorrectBias(IRF, w)

    % LP bias correction
    % "BCC" estimator in Herbst & Johanssen (2022)

    IRF_hor = length(IRF)-1;
    T = size(w,1);

    % ACF term in bias correction
    acf_corr = nan(1,IRF_hor);
    w = w - mean(w); % de-mean
    Sigma_0 = cov(w); % var-cov
    for j=1:IRF_hor
        Sigma_j = (w(1:end-j,:)'*w(j+1:end,:))/(T-j-1); % j-th autocov
        acf_corr(j) = 1+trace(Sigma_0\Sigma_j);
    end

    % Iterate on bias correction
    IRF_corr = IRF;
    for h=1:IRF_hor
        IRF_corr(h+1) = IRF(h+1) + (1/(T-h))*acf_corr(1:h)*IRF_corr(h:-1:1)';
    end

end