function obj = obj_fn_aux(poly,the_bias,the_std)

% get polynomial coefficients

a = poly(1);
b = poly(2);
c = poly(3);

% get predicted values

the_std_pred = a + b * the_bias + c * the_bias.^2;

% get error

obj = sum((the_std_pred - the_std).^2);