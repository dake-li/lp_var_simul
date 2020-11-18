function [MSE, BIAS2, VCE] = irf_perform_summary(est_irf,target_irf,settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    
    MSE.(thisMethod) = squeeze(mean((est_irf.(thisMethod) ...
        - permute(target_irf,[1 3 2])).^2, 2));
    
    BIAS2.(thisMethod) = (squeeze(mean(est_irf.(thisMethod) ...
        , 2)) - target_irf).^2;
    
    VCE.(thisMethod) = squeeze(var(est_irf.(thisMethod), 0, 2));
    
end
end

