function [data_sim_select] = select_data_fn(data_sim_all,settings,i_spec)
% Function for extracting the simulated data of \bar{w}_t from the
% simulated data of y_t, in a specific DGP

var_select = settings.specifications.var_select;
n_spec = settings.specifications.n_spec;

with_IV = settings.est.with_IV;

if with_IV == 1
    rho_select_grid_idx = settings.specifications.rho_select_grid_idx;
    this_rho_grid_idx = rho_select_grid_idx(i_spec);
else
    this_rho_grid_idx = 1;
end

data_sim_select.data_y = data_sim_all.data_y(:,var_select(i_spec,:));
data_sim_select.data_z = data_sim_all.data_z(:,this_rho_grid_idx);
data_sim_select.data_shock = data_sim_all.data_shock;

end

