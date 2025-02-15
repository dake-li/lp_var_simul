function IRF = compute_irfs(ABCD,shock_weight,IRF_hor);
% Function for computing the true IRF
    % Use a general ABCD representation of the encompassing model (DFM, DSGE or others):
        % state transition:  s_t = A * s_{t-1} + B * e_t
        % measurement eq:    y_t = C * s_{t-1} + D * e_t

    % Computes IRF with respect to the shock (shock_weight'*e_t)

% unpack settings

A = ABCD.A;
B = ABCD.B;
C = ABCD.C;
D = ABCD.D;
n_y = size(C,1);

% compute IRFs for all structural shocks

IRF = NaN(IRF_hor, n_y);
IRF(1,:) = (D * shock_weight)';

% go through horizon 1 to IRF_hor - 1
for i = 2:IRF_hor
    IRF(i,:) = C * A^(i-2) * B * shock_weight; % iterate to get IRF
end

end