function [B, Xb, W, Y_resw, Xb_resw] = locproj_design(y, x, w, H_min, H_max)
% Function for designing data matrix in penalized LP

    %%% Basic %%%
    % T:  sample size
    % K:  number of B-spline basis functions
    % HR: number of horizons
    % p:  number of controls
    
    %%% Input %%%
    % y:       response variable
    % x:       impulse variable
    % w:       controls (contemperaneous and lagged)
    % H_min:   minimum horizon
    % H_max:   maximum horizon
    
    %%% Output %%%
    % B:       B-spline values for each basis function at each horizon
    % Xb:      x * B
    % W:       controls collected for each horizon
    % Y_resw:  y residualized by w
    % Xb_resw: Xb residualized by w
    
    % prepare
    T  = length(y);
    HR = H_max + 1 - H_min;

    % construct the B-spline basis functions
    B = bspline( (H_min:H_max)' , H_min , H_max+1 , H_max+1-H_min , 3 );
    K = size( B , 2 );

    % Residualize Y and Xb on the unpenalized controls
    w = [ ones(T,1) w ]; % Add intercept to controls
    p = size(w,2);

    % placeholder for output matrices
    Xb = nan(T,K,HR);
    W = nan(T,p,HR);
    Y_resw = nan(T,HR);
    Xb_resw = nan(T,K,HR);

    % go thru each horizon
    for ih=1:HR

        h = H_min-1+ih;

        % Shift x and w relative to y, according to horizon h
        the_lagm = lagmatrix([x w], h);
        Xb(:,:,ih) = the_lagm(:,1)*B(ih,:); % x*B
        W(:,:,ih) = the_lagm(:,2:end); % w
        
        if nargout>3

            % Residualize y and x*B(h) on controls
            the_select = ~isnan(the_lagm(:,1));
            the_reg = the_lagm(the_select,2:end);
            the_outc = [y(the_select) Xb(the_select,:,ih)];
            the_beta = the_reg\the_outc;
            the_resw = the_outc - the_reg*the_beta;

            Y_resw(the_select,ih) = the_resw(:,1);
            Xb_resw(the_select,:,ih) = the_resw(:,2:end);
        
        end

    end

end


function B = bspline(x, xl, xr, ndx, bdeg)
    dx = (xr - xl) / ndx;
    t = xl + dx * [-bdeg:ndx-1];
    T = (0 * x + 1) * t;
    X = x * (0 * t + 1);
    P = (X - T) / dx;
    B = (T <= X) & (X < (T + dx));
    r = [2:length(t) 1];
    for k = 1:bdeg
        B = (P .* B + (k + 1 - P) .* B(:, r)) / k;
    end
end
