function res_new = compute_cum_irf(res, cum_irf_by_trans_code)
% Function to compute cumulative IRF for certain variables

% unpack
spec_id = res.spec_id;
dgp_type = res.dgp_type;
estimand_type = res.estimand_type;
lag_type = res.lag_type;
DFM_estimate = res.DFM_estimate;
DF_model = res.DF_model;
settings = res.settings;
results = res.results;

% decide which specfication needs cumulative IRF
resp_var_select = settings.specifications.var_select(:, settings.est.IRF_response_var_pos);
trans_code_for_resp_by_spec = DF_model.trans_code(resp_var_select);
cum_irf_by_spec = logical(cum_irf_by_trans_code(trans_code_for_resp_by_spec))';
cum_irf_across_all_var = logical(cum_irf_by_trans_code(DF_model.trans_code))';

% compute cumulative numbers for true IRF
DF_model.irf(:, cum_irf_across_all_var) = cumsum(DF_model.irf(:, cum_irf_across_all_var), 1);
DF_model.target_irf(:, cum_irf_by_spec) = cumsum(DF_model.target_irf(:, cum_irf_by_spec), 1);
if strcmp(estimand_type, 'Recursive')
    DF_model.VAR_irf(:, cum_irf_by_spec) = cumsum(DF_model.VAR_irf(:, cum_irf_by_spec), 1);
else
    DF_model.normalized_irf(:, cum_irf_by_spec) = cumsum(DF_model.normalized_irf(:, cum_irf_by_spec), 1);
end

% compute cumulative numbers for estimated IRF
method_names = fieldnames(results.irf);
for i = 1:length(method_names)
    results.irf.(method_names{i})(:,:, cum_irf_by_spec) = cumsum(results.irf.(method_names{i})(:,:, cum_irf_by_spec), 1);
end

% update settings
settings.est.cum_irf_by_trans_code = cum_irf_by_trans_code;
settings.specifications.cum_irf_by_spec = cum_irf_by_spec;

% pack up
res_new.spec_id = spec_id;
res_new.dgp_type = dgp_type;
res_new.estimand_type = estimand_type;
res_new.lag_type = lag_type;
res_new.DFM_estimate = DFM_estimate;
res_new.DF_model = DF_model;
res_new.settings = settings;
res_new.results = results;

end

