function ABCD = ABCD_fun_DFM(model)
% Function for transforming the encompassing DFM model into a general ABCD representation
    % Encompassing DFM model:
        % factor transition: f_t = \Phi(L) f_{t-1} + H \epsilon_t
        % observables:       X_t = \Lambda f_t + v_t
        % measurement error: v_{it} = \Delta_i(L) v_{i,t-1} + \Xi_i \xi_{it}

        % where f_t are factors, X_t are the full set of observables
        %       \epsilon_t are the structural shocks
        %       \xi_{it} are innovations for idiosyncratic measurement errors
    
    % ABCD representation:
        % state transition:  s_t = A * s_{t-1} + B * e_t
        % measurement eq:    y_t = C * s_{t-1} + D * e_t

        % where s_t = (f_t',...,f_{t-p_f+1}',v_t',...,v_{t-p_v+1}')', e_t = (epsilon_t',xi_t')'. 
        %       \epsilon_t are the structural shocks.
        %       \xi_t are innovations in measurement errors.
        %       y_t are observables, s_t are states. Warning: y_t correspond to X_t in our paper

    % Transforming formula can be found in Technical Companion Note
    
    % Input:
        % model: struct that contains all parameters in DFM
    % Output:
        % ABCD: struct that contains all the matrices of A B C D

aux = zeros(model.n_y,model.n_y*model.n_lags_uar);
for l=1:model.n_lags_uar
    aux(:,model.n_y*(l-1)+1:model.n_y*l) = diag(model.delta(:,l));
end

H = chol(model.Sigma_eta, 'lower');

ABCD.A = blkdiag(model.Phi, ...
                [aux; ...
                eye(model.n_y*(model.n_lags_uar-1)) zeros(model.n_y*(model.n_lags_uar-1),model.n_y)]);
ABCD.B = [H zeros(model.n_fac,model.n_y);
          zeros(model.n_fac*(model.n_lags_fac-1),model.n_fac+model.n_y);
          zeros(model.n_y,model.n_fac) diag(model.sigma_v);
          zeros(model.n_y*(model.n_lags_uar-1),model.n_fac+model.n_y)];
ABCD.C = [model.Lambda*model.Phi(1:model.n_fac,:) aux];
ABCD.D = [model.Lambda*H diag(model.sigma_v)];

end

