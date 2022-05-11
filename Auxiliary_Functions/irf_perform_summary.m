function [MSE, BIAS2, VCE] = irf_perform_summary(est_irf,target_irf,settings)
% Function for summarizing IRF estimates across simulations, in terms of MSE,
% bias squared and variance

% go thru each estimator
for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    
    BIAS2.(thisMethod) = (squeeze(est_irf.(thisMethod)(:,1,:)) - target_irf).^2;
    
    VCE.(thisMethod) = squeeze(est_irf.(thisMethod)(:,2,:)).^2;
    
    MSE.(thisMethod) = BIAS2.(thisMethod) + VCE.(thisMethod);
    
end

end

