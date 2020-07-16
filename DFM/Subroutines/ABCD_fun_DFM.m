function ABCD = ABCD_fun_DFM(model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ABCD.A = model.Phi;
ABCD.B = [chol(model.Sigma_eta, 'lower'); zeros(model.n_fac * (model.n_lags_fac - 1), model.n_fac)];
ABCD.C = kron([1, zeros(1, model.n_lags_fac - 1)], model.Lambda) * ABCD.A;
ABCD.D = kron([1, zeros(1, model.n_lags_fac - 1)], model.Lambda) * ABCD.B;
ABCD.E = zeros(model.n_w * model.n_lags_uar);
for ilag = 1:model.n_lags_uar
    ABCD.E(1:model.n_w, (ilag - 1) * model.n_w + (1:model.n_w)) = diag(model.delta(:, ilag));
end
ABCD.E((model.n_w + 1):end, 1:((model.n_lags_uar - 1) * model.n_w)) = eye((model.n_lags_uar - 1) * model.n_w);
ABCD.F = [diag(model.sigma_v); zeros((model.n_lags_uar - 1) * model.n_w, model.n_w)];

end

