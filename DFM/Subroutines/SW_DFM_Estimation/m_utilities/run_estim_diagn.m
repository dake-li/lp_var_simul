% DFM estimation diagnostics

quants = [0 0.1 0.25 0.5 0.75 0.9 1]; % Cross-sectional quantiles to display
maxlag = 10; % Max lag for info criteria

ns = length(DFM_estimate.sigma_v); % No. of series


%% Factors

[BIC,AIC] = IC_VAR(DFM_estimate.fac,maxlag);
figure;
plot([BIC' AIC']);
legend({'BIC','AIC'},'Location','NorthWest');
title('Factors: lag length info criteria');

disp('Factors: VAR eigenvalues');
disp(sort(abs(eig(DFM_estimate.Phi)),'descend')');

disp('Factors: innovation var-cov eigenvalues');
disp(sort(abs(eig(DFM_estimate.Sigma_eta)),'descend')');

disp('Factors: std dev');
disp(std(DFM_estimate.fac,'omitnan'));

disp('Factors: quantiles of R^2, in levels');
disp(quantile(DFM_estimate.r2,quants));

% Compute R^2 on differenced data
dat = DFM_estimate.fac*DFM_estimate.Lambda'+DFM_estimate.v'; % Reconstruct data
dat_diff = diff(dat,1);
fac_diff = diff(DFM_estimate.fac,1);
r2_diff = nan(ns,1);
for j=1:ns
    the_smpl = ~isnan(dat_diff(:,j));
    the_load_diff = fac_diff(the_smpl,:)\dat_diff(the_smpl,j);
    the_res_diff = dat_diff(the_smpl,j)-fac_diff(the_smpl,:)*the_load_diff;
    r2_diff(j) = 1-var(the_res_diff,'omitnan')/var(dat_diff(:,j),'omitnan');
end
disp('Factors: quantiles of R^2, in changes');
disp(quantile(r2_diff,quants));


%% Loadings

disp('Loadings: quantiles, by factor');
disp(quantile(DFM_estimate.Lambda,quants,1)');

disp('Loadings: absolute, quantiles, by factor');
disp(quantile(abs(DFM_estimate.Lambda),quants,1)');


%% Idiosyncratic residuals

% Compute summaries for each series
idio_roots = nan(ns,size(DFM_estimate.delta,2));
idio_bic = nan(ns,maxlag);
idio_aic = nan(ns,maxlag);
for j=1:ns
    idio_roots(j,:) = roots([-DFM_estimate.delta(j,end:-1:1) 1]); % AR roots
    the_v = DFM_estimate.v(j,:);
    [idio_bic(j,:),idio_aic(j,:)] = IC_VAR(the_v(~isnan(the_v))', maxlag); % Info criteria
end
[~,idio_bic_p] = min(idio_bic,[],2);
[~,idio_aic_p] = min(idio_aic,[],2);

disp('Idiosyncratic: quantiles of lag lengths selected, BIC');
disp(quantile(idio_bic_p,quants));

disp('Idiosyncratic: quantiles of lag lengths selected, AIC');
disp(quantile(idio_aic_p,quants));

disp('Idiosyncratic: quantiles of largest AR root');
disp(quantile(max(1./abs(idio_roots),[],2),quants));

disp('Idiosyncratic: quantiles of AR innovation std dev');
disp(quantile(DFM_estimate.sigma_v,quants));
