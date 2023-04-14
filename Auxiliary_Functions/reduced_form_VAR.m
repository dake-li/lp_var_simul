function out = reduced_form_VAR(ABCD,num_lag)

    % Compute reduced-form VAR representation
    % implied by ABCD model of the form
        % s_t = A*s_{t-1} + B*e_t
        % y_t = C*s_{t-1} + D*e_t
    % with e_t ~ WN(0,I)

    % First we compute the innovation representation
        % x_t = A*x_{t-1} + K*u_t
        % y_t = C*x_{t-1} + u_t
    % with x_t = E[y_t | y_{t-1},y_{t-2},...] and u_t ~ WN(0,Sigma),
    % as in Fernandez-Villaverde, Rubio-Ramirez, Sargent & Watson (AER 2007)

    % Then we compute the VAR(infinity) representation

    A = ABCD.A;
    B = ABCD.B;
    C = ABCD.C;
    D = ABCD.D;

    % Compute steady-state conditional variance in Kalman filter

    out.cond_var = cond_var_state(A, B, C, D);

    % Compute var-cov matrix of reduced-form innovations

    out.innov_var = C * out.cond_var * C' + D * D';

    % Compute Kalman gain

    K = (A * out.cond_var * C' + B * D') / out.innov_var;

    % Store innovation representation coefficients as ABCD model

    out.innov_ABCD.A = A;
    out.innov_ABCD.B = K;
    out.innov_ABCD.C = C;
    out.innov_ABCD.D = eye(size(C,1));

    % Compute coefficients in VAR(infinity) representation
    % See Hansen & Sargent (2014 book), chapter 8.6

    out.coef = cell(1,num_lag);
    AmKC = A-K*C;
    aux = C;
    for i=1:num_lag
        out.coef{i} = aux*K; % VAR coefs at lag i
        aux = aux*AmKC;
    end

end