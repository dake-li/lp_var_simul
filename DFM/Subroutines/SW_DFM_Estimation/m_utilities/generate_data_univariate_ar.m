function [y_mat] = generate_data_univariate_ar(T,nrep,ar_coef,se_ar,n_initial)
% 

% Parameters
nar = max(size(ar_coef));

% Initial values set to zero
y_mat = zeros(T+n_initial,nrep);
for t = nar+1:T+n_initial;
    y_mat(t,:) = se_ar*randn(1,nrep);
    for iar = 1:nar;
        y_mat(t,:) = y_mat(t,:)+ar_coef(iar)*y_mat(t-iar,:);
    end;
end;
y_mat = y_mat(n_initial+1:end,:);

    
end

