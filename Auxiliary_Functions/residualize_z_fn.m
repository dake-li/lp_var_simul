function [data_sim_new] = residualize_z_fn(data_sim,settings)
% Function for residualizing IV using lagged IV and endogenous variables w_{t-l}
    % This function should be used before running SVAR-IV estimators

% preparations

run('Estimation_Setup'); % common setup for all estimation methods

% estimate VAR

[~,~,~,~,Y_Res] = VAR(Y,nlags);

Z_Res = Y_Res(:,1);

data_sim_new.data_z = Z_Res;

data_sim_new.data_y = data_sim.data_y((nlags+1):end, :);
data_sim_new.data_shock = data_sim.data_shock((nlags+1):end);

end

