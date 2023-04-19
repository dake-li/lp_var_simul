function [dcov, dlrv] = dcov_dlrv(ABCD,settings)

    % Var-cov matrix and long-run var-cov matrix of dy_t = y_t - y_{t-1} in the model
        % s_t = A*s_{t-1} + B*e_t
        % y_t = C*s_{t-1} + D*e_t
    % with e_t ~ WN(0,I)

    A = ABCD.A;
    B = ABCD.B;
    C = ABCD.C;
    D = ABCD.D;
    n_s = size(A,1);
    lag_hor = settings.est.VAR_infinity_truncate;
    
    % Compute Cov(dy_t)

    dcov = D * D' + (C * B - D) * (C * B - D)';
    for i = 0:lag_hor
        aux = C * (A - eye(n_s)) * A^i * B;
        dcov = dcov + aux * aux';
    end

    % Compute LRV(dy_t)
    aux = D + (C * B - D);
    for i = 0:lag_hor
        aux = aux + C * (A - eye(n_s)) * A^i * B;
    end
    dlrv = aux * aux';

end