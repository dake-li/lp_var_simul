function [irf_draws] = IRF_BVAR(VAR_coef_post_mean,VAR_coef_post_vce_inv,ShockVector,nhorizons,ndraw)
% Auxiliary function for do posterior sampling of IRFs

nv = size(ShockVector,1);
nx = size(VAR_coef_post_mean,1) / nv;
nlags = (nx - 1) / nv;
irf_draws = NaN(nv,nhorizons + 1,ndraw);
VAR_coef_post_vce_sqrt = chol(inv(VAR_coef_post_vce_inv), 'lower');

% run multiple posterior draws
for idraw = 1:ndraw
    
    % make draws on VAR coef.
    this_VAR_coef = VAR_coef_post_mean + VAR_coef_post_vce_sqrt * randn(nx*nv, 1);
    
    % reshape to get By (reduced-form coef. on endogenous variables)
    this_Beta = reshape(this_VAR_coef, [nx, nv]);
    this_By = reshape(this_Beta(2:end,:),[nv,nlags,nv]); % lagged term
    this_By = permute(this_By,[3,1,2]);
    
    % compute IRF based on this draw
    irf_draws(:,:,idraw) = IRF_SVAR(this_By,ShockVector,nhorizons);
    
end

end

