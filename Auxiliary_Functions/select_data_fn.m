function [data_sim_select] = select_data_fn(data_sim_all,settings,i_spec)
% Function for extracting the simulated data of \bar{w}_t from the
% simulated data of y_t, in a specific DGP

var_select = settings.specifications.var_select;
n_spec = settings.specifications.n_spec;

with_IV = settings.est.with_IV;

if with_IV == 1
    IV_persistence_scale = settings.est.IV.IV_persistence_scale;
    n_rho_grid = length(IV_persistence_scale);
else
    n_rho_grid = 1;
end

n_spec_unique = n_spec / n_rho_grid;
i_rho_grid = ceil(i_spec / n_spec_unique);

data_sim_select.data_y = data_sim_all.data_y(:,var_select(i_spec,:));
data_sim_select.data_z = data_sim_all.data_z(:,i_rho_grid);
data_sim_select.data_shock = data_sim_all.data_shock;

end

