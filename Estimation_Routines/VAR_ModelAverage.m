function [combination_irf,weights,H] = VAR_ModelAverage(dat,recurShock,respV,p,hmax,options)
% Creates h-step impulse response functions from estimated VAR(p) model with intercept
% (modified based on Bruce Hansen's "cvar_ir.m" code)
% https://www.ssc.wisc.edu/~bhansen/progs/var.html

% Inputs:
%	dat 			nxm data matrix
%	p			VAR order p >= 1
%	hmax			horizon, hmax >= 1
%   recurShock  which shock?
%   respV       which respones?
%   options     numerical options for quadprog

% Outputs:
%	combination_irf		(hmax+1) orthogonalized impulse responses, horizon 0 to hmax
%	weights			M x hmax 	model weights, by horizon (1 to hmax) and response variable
%   H               m x m       choleskey decomp. of residual cov-var matrix


%% Prepare

% Dimensions
n = size(dat,1)-p;
m = size(dat,2);
k = m*p+1;

% Submodels
M = 2*p;
submodel_irf = zeros(hmax+1, M); 
combination_irf = zeros(hmax+1, 1); 


%% Full model estimation

% Least squares
[~,By,sigma,Q,e,Bt,y,x] = VAR(dat,p);
sel = [2:k 1]; % Swap constant and slope parameters
Bt = Bt(sel,:);
x = x(:,sel);
Q = Q(sel,sel);
B = Bt';
theta = Bt(:);
H = chol(sigma,'lower');

% Asymptotic variance
xe = repmat(x,1,m) .* e(:,kron(1:m,ones(1,k)));
omega = (xe'*xe) / (n-k);
WI = kron(eye(m),inv(Q));
V = WI*omega*WI;

% IRF
irf = IRF_SVAR(By,H(:,recurShock),hmax);
var_ir = irf(respV,2:end);
submodel_irf(:,M) = irf(respV,:);
combination_irf(1,1) = irf(respV,1);

% Jacobian
J = [eye(m*(p-1)), zeros(m*(p-1),m+1); zeros(1,m*p), 1];
P = [B(:,1:m*p), zeros(m,1); J];
Ph0 = [H; zeros(k-m,m)];
G0 = zeros(k*m,hmax);
for h = 1:hmax
    Gh = zeros(k*m,m^2);
    T2 = Ph0;
    for i = 1:h
      T1 = (P^(h-i))';
      Gh = Gh + kron(T1(1:m,1:m),T2);
      T2 = P*T2;
    end
    G0(:,h) = Gh(:, (respV - 1) * m + recurShock); % keep only response of interest
end


%% Estimate submodels and their loss

% Auxiliary matrices in loss function
V_G0 = V*G0;
WI_G0 = reshape(Q\reshape(G0, [k m*hmax]), [k*m hmax]); % Acklam section 10.1.8

% Loop over submodels
ir_diff = zeros(M,hmax);
Kr = zeros(M,hmax);

for r = 1:M-1
    
    % Identify constrained parameters and estimate constrained model
    B_sel = [false(m,m*p) true(m,1)]; % Default: only intercepts are unconstrained
    if r<=p % AR(1) to AR(p)
        B_sel(:,1:m*r) = repmat(eye(m)==1,1,r);
    else % VAR(1) to VAR(p-1)
        r1 = r-p;
        B_sel(:,1:m*r1) = true;
    end
    B_sel_t = B_sel';
    sel_vec = ~B_sel_t(:); % Constrained indices in theta=vec(B')
    WI_sel = WI(sel_vec,sel_vec);
    
    thetar = theta - WI(:,sel_vec)*(WI_sel\theta(sel_vec));
    Br = (reshape(thetar,k,m))';
    er = y - x*Br';
    sigmar = (er'*er) / (n-k+sum(sel_vec)/m);
    Hr = chol(sigmar,'lower');
    Pr = [Br(:,1:(m*p)), zeros(m,1); J];
    Ph = [Hr; zeros(k-m,m)];
    submodel_irf(1,r) = Hr(respV, recurShock);
    
    % Loss
    WI_sel_inv_V_G0 = WI_sel\V_G0(sel_vec,:);
    for h = 1:hmax
        Ph = Pr*Ph;
        Phm = Ph(respV, recurShock);
        diff = Phm - var_ir(1,h);
        ir_diff(r,h) = diff;
        submodel_irf(h+1,r) = Phm;
        Kr(r,h) = WI_G0(sel_vec,h)'*WI_sel_inv_V_G0(:,h);
    end
    
end


%% Compute weights for each horizon

weights = zeros(M, hmax);
ub = ones(M,1);
lb = zeros(M,1);
for h = 1:hmax
    ird = ir_diff(:,h);
    J0 = n*(ird*ird');
    K = - Kr(:,h)';
    w = quadprog(J0,K,[],[],ub',1,lb,ub,[],options);
    weights(:,h) = w;
    combination_irf(h+1,1) = submodel_irf(h+1,:)*w;
end

end

