function std_frontier_grid = std_frontier_fn(bias,std,bias_grid)

    % Solve for "largest" quadratic function that lies below all points
    
    % First remove any dominated points, as we don't want those to
    % influence the lower envelope
    bias_diff = bias-bias'; % Pairwise differences of bias
    std_diff = std-std'; % Pairwise differences of stdev
    dom = all(bias_diff>=0 & std_diff>=0, 2); % Dominated procedures
    bias = bias(~dom);
    std = std(~dom);
    
    % Solve for quadratic lower "envelope"
    X = bias.^(0:2);
    coefs = linprog(-sum(X)', ... % Maximize average value of quadratic function on grid points
                    X, std, ... % ... subject to the function lying below all data points
                    [], [], ...
                    [-Inf -Inf 0], ... % Enforce that the quadratic polynomial is "smiling"
                    []);
    
    % Evaluate quadratic function on grid
    X_grid = bias_grid.^(0:2);
    std_frontier_grid = X_grid*coefs;


end