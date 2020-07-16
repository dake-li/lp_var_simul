function [IRF, shock_weight] = compute_irfs(model,settings);

% prepare

IRF_hor = settings.est.IRF_hor;
manual_shock_pos = settings.est.manual_shock_pos;
estimate_shock_weight = settings.est.estimate_shock_weight;
shock_weight_calibrate = settings.est.shock_weight_calibrate;
shock_optimize_var_IRF = settings.est.shock_optimize_var_IRF;
n_y = model.n_y;
n_eps = model.n_eps;

A = model.ABCD.A;
B = model.ABCD.B;
C = model.ABCD.C;
D = model.ABCD.D;

if (estimate_shock_weight==1) && (shock_weight_calibrate==1)
    calibrated_shock_weight = model.calibrated_shock_weight;
end

% compute IRFs for all structural shocks

IRF = NaN(IRF_hor, n_y);

for i = 1:IRF_hor
    if i == 1
        resp = D;
        if (estimate_shock_weight==1) && (shock_weight_calibrate==0)
            shock_weight = resp(shock_optimize_var_IRF,:)'; % optimal weight in the same direction as impulse
            shock_weight = shock_weight / sqrt(sum(shock_weight.^2)); % normalize weight
        elseif (estimate_shock_weight==1) && (shock_weight_calibrate==1)
            shock_weight = calibrated_shock_weight; % calibrated shock weights
        else
            shock_weight = zeros(n_eps, 1); % manually choose shock
            shock_weight(manual_shock_pos) = 1;
        end
        IRF(i,:) = (resp * shock_weight)';
    else
        IRF(i,:) = C * A^(i-2) * B * shock_weight;
    end
end

end