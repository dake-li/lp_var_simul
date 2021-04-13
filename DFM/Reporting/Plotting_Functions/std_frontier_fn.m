function std_frontier = std_frontier_fn(the_bias,the_std,bias_frontier);

obj_fn    = @(poly) obj_fn_aux(poly,the_bias,the_std);
constr_fn = @(poly) constr_fn_aux(poly,the_bias,the_std,bias_frontier);

solution = fmincon(obj_fn,[0 0 0],[],[],[],[],[],[],constr_fn);
% solution = fmincon(obj_fn,[0 0 0],[],[],[],[],[],[],[]);

std_frontier  = solution(1) + solution(2) * bias_frontier + solution(3) * bias_frontier.^2;