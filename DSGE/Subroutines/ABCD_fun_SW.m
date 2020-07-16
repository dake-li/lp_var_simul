function ABCD = ABCD_fun_SW(model);

% prepare inputs

decision = model.decision;
n_y      = model.n_y;
n_s      = model.n_s;
n_eps    = model.n_eps;
obs      = model.obs;

% ABCD representation in macro aggregates y

A = decision(1:n_s,1:n_s); A = A';
B = decision(n_s+1:n_s+n_eps,1:n_s); B = B';
C = decision(1:n_s,obs); C = C';
D = decision(n_s+1:n_s+n_eps,obs); D = D';

% collect results

ABCD.A = A;
ABCD.B = B;
ABCD.C = C;
ABCD.D = D;
ABCD.E = zeros(n_y);
ABCD.F = zeros(n_y);

end