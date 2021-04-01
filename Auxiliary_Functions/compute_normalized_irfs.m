function IRF = compute_normalized_irfs(model,settings);
% Function for computing the true normalized IRF in the ObsShock and IV experiments

% unpack
true_irf = model.irf;
var_select = settings.specifications.var_select;
response_pos = settings.est.IRF_response_var_pos;
normalize_pos = settings.est.est_normalize_var_pos;
with_shock = settings.est.with_shock;
recursive_shock = settings.est.recursive_shock;
with_IV = settings.est.with_IV;
if with_shock == 1
    normalize_with_shock_std_dev = settings.est.normalize_with_shock_std_dev;
end

% pick IRF
IRF_response = true_irf(:, var_select(:, response_pos));

% normalize with the unit-std-dev of shock or the unit of normalization variable
if (with_shock == 1) && (normalize_with_shock_std_dev == 1)
    IRF = IRF_response;
else
    IRF = IRF_response ./ true_irf(1, var_select(:, normalize_pos));
end

end

