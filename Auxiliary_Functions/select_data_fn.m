function [data_sim_select] = select_data_fn(data_sim_all,settings,i_spec)
% Function for extracting the simulated data of \bar{w}_t from the
% simulated data of y_t, in a specific DGP

var_select = settings.specifications.var_select;

data_sim_select.data_y = data_sim_all.data_y(:,var_select(i_spec,:));
data_sim_select.data_z = data_sim_all.data_z;
data_sim_select.data_shock = data_sim_all.data_shock;

end

