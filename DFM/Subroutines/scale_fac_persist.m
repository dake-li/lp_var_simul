function DF_model_new = scale_fac_persist(DF_model)
% Function for scaling up factor persistence in the encompassing DFM

% unpack

n_fac = DF_model.n_fac;
n_lags_fac = DF_model.n_lags_fac;
Phi = DF_model.Phi;
Sigma_eta = DF_model.Sigma_eta;
target_largest_root = DF_model.fac_persist.target_largest_root;
keep_fac_variation = DF_model.fac_persist.keep_fac_variation;

% prepare scaling up autocorrelation coefficient Phi

DF_model_new = DF_model;
Phi_new = Phi;
old_largest_root = max(abs(eig(Phi)));
scale_Phi = target_largest_root / old_largest_root;

% scale up Phi

for i_lag = 1:n_lags_fac
    Phi_new(1:n_fac, (i_lag-1)*n_fac + (1:n_fac)) = Phi_new(1:n_fac, (i_lag-1)*n_fac + (1:n_fac)) * (scale_Phi^i_lag);
end

% check if keep factor variation the same

if keep_fac_variation == 1

    % use companion form to compute var-cov of factors with old and new Phi

    Sigma_eta_comp = zeros(n_fac*n_lags_fac);
    Sigma_eta_comp(1:n_fac, 1:n_fac) = Sigma_eta;
    
    var_cov_fac_comp = (eye((n_fac*n_lags_fac)^2) - kron(Phi,Phi)) \ Sigma_eta_comp(:);
    var_cov_fac_comp = reshape(var_cov_fac_comp, [n_fac*n_lags_fac, n_fac*n_lags_fac]);
    var_cov_fac = var_cov_fac_comp(1:n_fac, 1:n_fac);

    var_cov_fac_comp_new = (eye((n_fac*n_lags_fac)^2) - kron(Phi_new,Phi_new)) \ Sigma_eta_comp(:);
    var_cov_fac_comp_new = reshape(var_cov_fac_comp_new, [n_fac*n_lags_fac, n_fac*n_lags_fac]);
    var_cov_fac_new = var_cov_fac_comp_new(1:n_fac, 1:n_fac);

    % prepare scaling down factor innovation Sigma_eta

    scale_Sigma = trace(var_cov_fac) / trace(var_cov_fac_new);

    % scale down Sigma_eta

    Sigma_eta_new = Sigma_eta * scale_Sigma;

else

    % not adjust Sigma_eta

    Sigma_eta_new = Sigma_eta;

end

% pack

DF_model_new.Phi = Phi_new;
DF_model_new.Sigma_eta = Sigma_eta_new;

end