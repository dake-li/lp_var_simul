function obj = locproj(varargin)
% locproj
%
%   locporj(y,x,w,H_min,H_max,type)
%   locporj(y,x,w,H_min,H_max,type,r,lambda)
%
%   
    switch length(varargin)
        case 6
            y     = varargin{1};
            x     = varargin{2};
            w     = varargin{3};
            H_min = varargin{4};
            H_max = varargin{5};
            type  = varargin{6};

        case 8
            y       = varargin{1};
            x       = varargin{2};
            w       = varargin{3};
            H_min   = varargin{4};
            H_max   = varargin{5};
            type    = varargin{6};            
            r       = varargin{7};
            lambda  = varargin{8};
            
        otherwise
            error('wrong number of input arguments')
    end    

    obj = struct();
    
    if isempty(w)
        delta = std(x);
    else
        delta = std( x-w*inv(w'*w)*w'*x );
    end
    delta = 1; % interested in y's response when x increases by one
    
    isreg = strcmp('reg',type);
    
    T  = length(y);
    HR = H_max + 1 - H_min;
    
    % construct the B-spline basis functions
    if ~isreg
        B = bspline( (H_min:H_max)' , H_min , H_max+1 , H_max+1-H_min , 3 );
        K = size( B , 2 );
    else
        K = HR;
    end

    % building up the regression representation of the local projection
    idx = nan( (H_max+1)*T , 2 );
    Y   = nan( (H_max+1)*T , 1 );
    Xb  = zeros( (H_max+1)*T , K );
    Xc  = zeros( (H_max+1)*T , HR , size(w,2)+1 );
    % Xc  = zeros( (H_max+1)*T , K , size(w,2)+1 );
    
    w = [ ones(T,1) w ];
    
    for t = 1:T-H_min
        
        idx_beg = (t-1)*HR + 1;
        idx_end = t*HR;

        idx( idx_beg:idx_end , 1 ) = t;
        idx( idx_beg:idx_end , 2 ) = H_min:H_max;
        
        % y
        y_range = (t+H_min) : min((t+H_max),T)';
        Y( idx_beg:idx_end ) = [ y( y_range ) ; nan(HR-length(y_range),1) ];

        % x
        if isreg
            Xb( idx_beg:idx_end , : ) = eye(HR)*x(t);
        else
            Xb( idx_beg:idx_end , : ) = B*x(t);
        end

        % w
        for i = 1:size(w,2)
            Xc( idx_beg:idx_end , : , i ) = eye(HR)*w(t,i); % no basis fcn for control var
            % Xc( idx_beg:idx_end , : , i ) = B*w(t,i);
        end
        
    end
    
    X = Xb;
    for i = 1:size(w,2)
        X = [X Xc(:,:,i)];
    end
    
    select = isfinite(Y);  
    idx = idx(select,:);
    Y   = Y(select);
    X   = X(select,:);
    X   = sparse(X);

    % estimation
    IR  = zeros(H_max+1,1);
    
    if isreg

        theta     = ( X'*X )\( X'*Y );
        IR((H_min+1):end) = theta(1:K) * delta;
    
    else

        P = zeros( size(X,2) );

        D = eye(K);
        for k = 1:r 
            D = diff(D);
        end
        
        P(1:K,1:K) = D' * D;

        theta = ( X'*X + lambda*P )\( X'*Y );
        
        IR((1+H_min):end) = B * theta(1:K) * delta;
    end
    
    % pack everything up
    obj.T     = T;
    obj.H_min = H_min;
    obj.H_max = H_max;
    obj.HR    = HR;
    obj.K     = K;
    
    if isreg
        obj.B      = 0;
        obj.P      = zeros( size(X,2) );
        obj.lambda = 0;
    else
        obj.B      = B;
        obj.P      = P;
        obj.D      = D;
        obj.lambda = lambda;
    end
    
    obj.type  = type;
    obj.delta = delta;
    obj.idx   = idx;
    obj.Y     = Y;
    obj.X     = X;
    obj.theta = theta;
    obj.IR    = IR;
    obj.delta = delta;
    
    % debug stuff
    % Bs is for display / debugging pourposes only
    % Bs = bspline( (H_min:0.1:H)' , H_min , H+1 , H+1-H_min , 3 );
    % obj.Bs = Bs;
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
