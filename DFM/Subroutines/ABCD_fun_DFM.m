function ABCD = ABCD_fun_DFM(model)
% Function for transforming the encompassing DFM model into a general ABCDEF representation
    % Encompassing DFM model:
        % factor transition: f_t = \Phi(L) f_{t-1} + H \epsilon_t
        % observables:       X_t = \Lambda f_t + v_t
        % measurement error: v_{it} = \Delta_i(L) v_{i,t-1} + \Xi_i \xi_{it}

        % where f_t are factors, X_t are the full set of observables
        %       \epsilon_t are the structural shocks
        %       \xi_{it} are innovations for idiosyncratic measurement errors
    
    % ABCDEF representation:
        % state transition:  s_t = A * s_{t-1} + B * \epsilon_t
        % measurement eq:    y_t = C * s_{t-1} + D * \epsilon_t + e^*_t
        % measurement error: e_t = E * e_{t-1} + F * \omega_t

        % where e_t = (e^*_t', e^*_{t-1}', ...)'. Warning: e^*_t corresponds to v_t in our paper
        %       \epsilon_t are the structural shocks.
        %       \omega_t are innovations in measurement errors. Warning: \omega_t corresponds to \xi_t in our paper
        %       y_t are observables, s_t are states. Warning: y_t correspond to X_t in our paper

    % Transforming formula can be found in Technical Companion Note
    
    % Input:
        % model: struct that contains all parameters in DFM
    % Output:
        % ABCD: struct that contains all the matrices of A B C D E F
        
% linking DFM coefficients with ABCDEF
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

