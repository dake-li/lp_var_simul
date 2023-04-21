function dABCD = ABCD_diff(ABCD)

    % Derive ABCD representation of dy_t = y_t - y_{t-1} in the model
        % s_t = A*s_{t-1} + B*e_t
        % y_t = C*s_{t-1} + D*e_t
    % with e_t ~ WN(0,I)

    % The ABCD model for dy_t has the form
        % x_t  = A0*x_{t-1} + B0*e_t
        % dy_t = C0*x_{t-1} + D0*e_t
    % with x_t = (beta'*s_{t-1}, e_t) and A-I = alpha*beta'

    A = ABCD.A;
    B = ABCD.B;
    C = ABCD.C;
    D = ABCD.D;
    [n_s,n_e] = size(B);
    tol = 1e-4; % Numerical tolerance for unit roots
    
    % Reduced-rank decomposition of A-I=alpha*beta'
    [U,S,V] = svd(A-eye(n_s));
    nonzero_sv = (abs(diag(S))>tol);
    alpha = U(:,nonzero_sv)*S(nonzero_sv,nonzero_sv);
    beta = V(:,nonzero_sv);
    k = size(alpha,2);

    % Build ABCD model for dy_t
    dABCD.A = [[eye(k)+beta'*alpha, beta'*B]; zeros(n_e,k+n_e)];
    dABCD.B = [zeros(k,n_e); eye(n_e)];
    dABCD.C = [C*alpha, C*B-D];
    dABCD.D = D;
    
end