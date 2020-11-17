function data_sim = generate_data(model,settings)

% preparations

T      = settings.simul.T;
T_burn = settings.simul.T_burn;

A = model.ABCD.A;
B = model.ABCD.B;
C = model.ABCD.C;
D = model.ABCD.D;
E = model.ABCD.E;
F = model.ABCD.F;

n_s   = model.n_s;
n_eps = model.n_eps;
n_y   = model.n_y;
n_w   = model.n_w;
n_e   = model.n_e;

with_IV = settings.est.with_IV;

if with_IV == 1
    rho     = model.IV.rho;
    alpha   = model.IV.alpha;
    sigma_v = model.IV.sigma_v;
else % meaningless placeholders
    rho = 0.1;
    alpha = 1;
    sigma_v = 1;
end

shock_weight = settings.est.shock_weight;

% draw shocks

data_eps = randn(T_burn+T,n_eps);
data_w   = randn(T_burn+T,n_w);

% simulate states & measurement error

s = zeros(n_s,1);
e = zeros(n_e,1);
data_s = NaN(T_burn + T, n_s);
data_e = NaN(T_burn + T, n_e);
for t = 1:(T_burn + T)
    s = A * s + B * data_eps(t,:)';
    data_s(t,:) = s';
    e = E * e + F * data_w(t,:)';
    data_e(t,:) = e';
end

% simulate observables

data_y = NaN(T, n_y);
for t = 1:T
    data_y(t,:) = C * data_s(T_burn + (t-1),:)' + D * data_eps(T_burn + t,:)' + data_e(T_burn + t, 1:n_y)';
end

% simulate IV

z = filter(1, [1 -rho], alpha * data_eps * shock_weight + sigma_v * randn(T+T_burn,1));
data_z = z(T_burn+1:end);

% collect results and shift timing

data_sim.data_y     = data_y;
data_sim.data_shock = data_eps((T_burn+1):(T_burn+T), :) * shock_weight;
data_sim.data_z     = data_z;

end