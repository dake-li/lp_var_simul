function [IR, theta, gamma] = locproj(y, x, w, H_min, H_max, r, lambda)

    % Smooth/penalized local projection (Barnichon & Brownlees, 2019)
    
    % Compute design matrix
    [B, Xb, W, Y_resw, Xb_resw] = locproj_design(y, x, w, H_min, H_max);
    
    % Compute impulse responses and coefficients
    [IR, theta, gamma] = locproj_partitioned(y, B, Xb, W, Y_resw, Xb_resw, r, lambda);
    
end
