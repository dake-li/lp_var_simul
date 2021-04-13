function [constr_val,ceq] = constr_fn_aux(poly,the_bias,the_std,bias_frontier)

% get polynomial coefficients

a = poly(1);
b = poly(2);
c = poly(3);

% get predicted values

the_std_pred = a + b * the_bias + c * the_bias.^2;
std_frontier = a + b * bias_frontier + c * bias_frontier.^2;

% get error

constr_val = [the_std_pred - the_std;std_frontier(2:end) - std_frontier(1:end-1)]; % constraints: predicted values <= actual values, decreasing slope
ceq = [];