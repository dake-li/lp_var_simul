function IV_strength = IV_strength(ABCD, IV, sigma_v, shock_weight, VAR_infinity_truncate, normalize_var_pos)

% R^2 in regression of normalization variable i_t on IV z_t
% after controlling for lagged observables

% Augment ABCD representation with (residualized) IV tilde{z}_t = z_t-rho*z_{t-1}
[n_y,n_s] = size(ABCD.C);
ABCD_withIV.A = ABCD.A;
ABCD_withIV.B = [ABCD.B zeros(n_s,1)]; % Add IV measurement error in last column
ABCD_withIV.C = [ABCD.C; zeros(1,n_s)]; % Add IV in bottom row
ABCD_withIV.D = [ABCD.D zeros(n_y,1);
                 IV.alpha*shock_weight' sigma_v]; % Add IV in bottom row and measurement error in last column

% VAR(infinity) representation with IV
ABCD_withIV_small = ABCD_reduce(ABCD_withIV); % Reduce dimensionality
red_form_withIV = reduced_form_VAR(ABCD_withIV_small,VAR_infinity_truncate);
innov_var = red_form_withIV.innov_var; % Var-cov matrix of Wold innovations

% IV strength = squared correlation of Wold innovations for i_t and z_t
IV_strength = innov_var(normalize_var_pos,end)^2/(innov_var(normalize_var_pos,normalize_var_pos)*innov_var(end,end));

end