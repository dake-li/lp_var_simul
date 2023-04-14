function [cov, lrv] = cov_lrv(ABCD)

    % Var-cov matrix and long-run var-cov matrix of y_t in the model
        % s_t = A*s_{t-1} + B*e_t
        % y_t = C*s_{t-1} + D*e_t
    % with e_t ~ WN(0,I)

    A = ABCD.A;
    B = ABCD.B;
    C = ABCD.C;
    D = ABCD.D;
    n_s = size(A,1);
    
    % Compute Cov(y_t)
    BB = B*B';
    cov_s = reshape((eye(n_s^2)-kron(A,A))\BB(:),n_s,n_s);
    cov = C*cov_s*C' + D*D';

    % Compute LRV(y_t) using y_t = [C(I-AL)^{-1}BL+D]e_t.
    aux = C*((eye(n_s)-A)\B)+D;
    lrv = aux*aux';

end