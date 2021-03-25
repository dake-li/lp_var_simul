function IRF = compute_normalized_irfs(model,settings);
% Function for computing the true normalized IRF in the ObsShock and IV experiments

true_irf = model.irf;
var_select = settings.specifications.var_select;
response_pos = settings.est.IRF_response_var_pos;
normalize_pos = settings.est.est_normalize_var_pos;

IRF_response = true_irf(:, var_select(:, response_pos));
IRF = IRF_response ./ true_irf(1, var_select(:, normalize_pos));

end

