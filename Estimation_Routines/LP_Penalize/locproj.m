function [IR, theta, gamma] = locproj(y, x, w, H_min, H_max, r, lambda)
% Function for estimating IRF and coefficients in penalized LP
    % Reference: Smooth/penalized local projection (Barnichon & Brownlees, 2019)
    
    %%% Input %%%
    % y:       response variable
    % x:       impulse variable
    % w:       controls (contemperaneous and lagged)
    % H_min:   minimum horizon
    % H_max:   maximum horizon
    % r:       order of finite difference operator
    % lambda:  penalty strength
    
    %%% Output %%%
    % IR:    impulse response
    % theta: penalized coef, for Xb. Warning: correspond to b_k in our paper
    % gamma: unpenalized coef, for W. Warning: correspond to \zeta_h and \phi_{h,l} in our paper
    
    % Design data matrix
    [B, Xb, W, Y_resw, Xb_resw] = locproj_design(y, x, w, H_min, H_max);
    
    % Compute impulse responses and coefficients using partitioned formula
    [IR, theta, gamma] = locproj_partitioned(y, B, Xb, W, Y_resw, Xb_resw, r, lambda);
    
end
