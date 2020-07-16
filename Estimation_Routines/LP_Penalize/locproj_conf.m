function obj = locproj_conf(varargin)

    switch length(varargin)
        case 2
            obj     = varargin{1};
            H       = varargin{2};
            lambda  = 0.0;
        case 3
            obj     = varargin{1};
            H       = varargin{2};
            lambda  = varargin{3};
            
        otherwise
            error('wrong number of input arguments')
    end 

    % POINT ESTIMATE
    XXP   = ( obj.X'*obj.X + lambda * obj.P );
    theta = XXP \ ( obj.X'*obj.Y );

    % COMPUTE NW ESTIMATOR
    % BREAD
    bread = XXP^-1;

    % MEAT
    nlag    = H;
    T       = obj.T;
    npar    = length(obj.theta);
    weights = [ 0.5 (nlag+1-(1:nlag))/(nlag+1) ];    
    idx     = obj.idx;
    U       = obj.Y - obj.X * obj.theta;    
    X       = obj.X ;
    V       = zeros( npar , npar );

    for l = 0:nlag
        GplusGprime = zeros( npar , npar );
        for t = (l+1):(T-obj.HR-1) 
            S1 = X( idx(:,1)==t , : )' * U( idx(:,1)==t );
            S2 = X( idx(:,1)==(t-l) , : )' * U( idx(:,1)==(t-l) );
            GplusGprime = GplusGprime + S1 * S2' + S2 * S1';
        end
        V = V + weights(l+1) * GplusGprime;
    end
    meat = V;

    VC = bread * meat * bread;
    
    conf = nan( obj.H_max+1 , 2 );

    if strcmp(obj.type,'reg')==1
        ster = sqrt( diag( VC( 1:(obj.H_max+1-obj.H_min) , 1:(obj.H_max+1-obj.H_min) ) ) );                
        conf(1+obj.H_min:end,1) = theta(1:obj.K)*obj.delta + ster*obj.delta*norminv(0.05);
        conf(1+obj.H_min:end,2) = theta(1:obj.K)*obj.delta + ster*obj.delta*norminv(0.95);
    else
        ster = sqrt( diag( obj.B*VC( 1:obj.K , 1:obj.K )*obj.B' ) );                
        conf(1+obj.H_min:end,1) = obj.B*theta(1:obj.K)*obj.delta + ster*obj.delta*norminv(0.05);
        conf(1+obj.H_min:end,2) = obj.B*theta(1:obj.K)*obj.delta + ster*obj.delta*norminv(0.95);
    end

    obj.ster = ster;
    obj.conf = conf;

end
