function [IRF, shock_weight] = compute_irfs(model,settings);
% Function for defining the true shock and computing the true IRF
    % Use a general ABCDEF representation of the encompassing model (DFM, DSGE or others):
        % state transition:  s_t = A * s_{t-1} + B * \epsilon_t
        % measurement eq:    y_t = C * s_{t-1} + D * \epsilon_t + e^*_t
        % measurement error: e_t = E * e_{t-1} + F * \omega_t

        % where e_t = (e^*_t', e^*_{t-1}', ...)'. Warning: e^*_t corresponds to v_t in the paper
        %       \epsilon_t are the structural shocks.
        %       \omega_t are innovations in measurement errors. Warning: \omega_t corresponds to \xi_t in the paper
        %       y_t are observables. Warning: corresponds to X_t in the paper
    
    % In this function, we can compute the optimal (IRF-maximizing) linear
    % combination of these structural shocks (saved as shock_weight) and define it as the true shock
    % Note that shock_weight' * \epsilon_t is then equal to the true shock \epsilon_{1t}

% unpack settings

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

% check if shock weight has been estimated in IV calibration

if (estimate_shock_weight==1) && (shock_weight_calibrate==1)
    calibrated_shock_weight = model.calibrated_shock_weight;
end

% compute IRFs for all structural shocks

IRF = NaN(IRF_hor, n_y);

% go through horizon 0 to IRF_hor - 1
for i = 1:IRF_hor
    
    if i == 1
    
        resp = D;
        
        if (estimate_shock_weight==1) && (shock_weight_calibrate==0)
            
            % if want to estimate the optimal shock weight to maximize impact response
            shock_weight = resp(shock_optimize_var_IRF,:)'; % optimal weight is in the same direction as impulse response
            shock_weight = shock_weight / sqrt(sum(shock_weight.^2)); % normalize weight
            
        elseif (estimate_shock_weight==1) && (shock_weight_calibrate==1)
            
            % if want to use the shock weight computed in IV calibration
            shock_weight = calibrated_shock_weight; % calibrated shock weights
            
        else
            
            % if want to manually choose the true shock
            shock_weight = zeros(n_eps, 1); % manually choose shock
            shock_weight(manual_shock_pos) = 1;
            
        end
        
        % IRF to the true shock is now the linear combo of IRF to all the structural shocks
        IRF(i,:) = (resp * shock_weight)'; 
        
    else
        
        IRF(i,:) = C * A^(i-2) * B * shock_weight; % iterate to get IRF
        
    end
    
end

end