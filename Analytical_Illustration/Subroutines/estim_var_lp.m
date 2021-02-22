function [ir_estims,Ahat_2,kappahat] = estim_var_lp(Y,E,h)

    X = [E Y];
    X_lag = lagmatrix(X,1);
    ir_estims = nan(1,2);
    
    % VAR
    Ahat= (X_lag(2:end,:)\X(2:end,:))';
    res = X(2:end,:)-X_lag(2:end,:)*Ahat';
    Sigmahat = (res'*res)/(size(X,1)-1-size(X,2));
    Sigmahat_chol = chol(Sigmahat,'lower');
    gammahat = Sigmahat_chol(:,1)/Sigmahat_chol(1,1);
    ir_red = Ahat^h;
    ir_estims(1) = ir_red(2,:)*gammahat;
    Ahat_2 = Ahat(2,:);
    kappahat = gammahat(2);
    
    % LP
    Xt = [E X_lag];
    betahat = Xt(2:end-h,:)\Y(h+2:end);
    ir_estims(2) = betahat(1);

end