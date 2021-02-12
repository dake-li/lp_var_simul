function [IR, theta, gamma] = locproj_partitioned(y, B, Xb, W, Y_resw, Xb_resw, r, lambda)

    % Compute penalized local projection, exploiting partitioned formula

    % Uses results from  R. W. Farebrother (1978), Partitioned Ridge Regression, Technometrics, 20:2, 121-122
    % https://doi.org/10.1080/00401706.1978.10489635


    [T,K,HR] = size(Xb);
    p = size(W,2);
    
    % r-th difference matrix
    D = eye(K);
    for k = 1:r 
        D = diff(D);
    end

    % First compute penalized coefficients
    Xb_resw_stack = reshape(permute(Xb_resw, [1 3 2]), HR*T, K);
    Y_resw_stack = Y_resw(:);
    select = ~isnan(Y_resw_stack);

    theta = [Xb_resw_stack(select,:); sqrt(lambda)*D]\[Y_resw_stack(select); zeros(K-r,1)];

    IR = B * theta; % Impulse responses

    % Then compute unpenalized coefficients horizon by horizon
    gamma = nan(p,HR);
    if nargout>2
        for ih=1:HR
            the_select = ~isnan(Y_resw(:,ih));
            gamma(:,ih) = W(the_select,:,ih)\(y(the_select)-Xb(the_select,:,ih)*theta);
        end
    end
    
end
