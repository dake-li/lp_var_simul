function data_sim = generate_data(model,settings)
% Function for generating simulated data
    % Use a general ABCD representation of the encompassing model (DFM, DSGE or others):
        % state transition:  s_t = A * s_{t-1} + B * e_t
        % measurement eq:    y_t = C * s_{t-1} + D * e_t

    % Warning: shock_weight' * e_t corresponds to true shock \epsilon_{1t} in our paper
    
    % external IV: z_t = \rho * z_{t-1} + \alpha * (shock_weight' * \epsilon_t) + \nu_t

% unpack settings

T      = settings.simul.T;
T_burn = settings.simul.T_burn;

A = model.ABCD.A;
B = model.ABCD.B;
C = model.ABCD.C;
D = model.ABCD.D;

[n_s,n_e] = size(B);

with_IV = settings.est.with_IV;

if with_IV == 1
    rho_grid = model.IV.rho_grid;
    alpha    = model.IV.alpha;
    sigma_v  = model.IV.sigma_v;
else % meaningless placeholders in the case of no IV
    rho_grid = 0.1;
    alpha = 1;
    sigma_v = 1;
end

shock_weight = settings.est.shock_weight;

% draw shocks

data_e = randn(T_burn+T,n_e);

% simulate states & measurement error

s = zeros(n_s,1);
data_s = NaN(T_burn + T, n_s);
for t = 1:(T_burn + T)
    s = A * s + B * data_e(t,:)';
    data_s(t,:) = s';
end

% simulate observables

data_y = data_s(T_burn:end-1,:)*C' + data_e(T_burn+1:end,:)*D';

% simulate IV
z = NaN(T+T_burn, length(rho_grid)); % multiple IV persistence setups
for idx = 1:length(rho_grid)
    z(:,idx) = filter(1, [1 -rho_grid(idx)], alpha * data_e * shock_weight + sigma_v * randn(T+T_burn,1));
end

data_z = z(T_burn+1:end,:);

% collect results and shift timing

data_sim.data_y     = data_y;
data_sim.data_shock = data_e((T_burn+1):(T_burn+T), :) * shock_weight;
data_sim.data_z     = data_z;

end