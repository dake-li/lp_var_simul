function [combination_irf,weights,H] = VAR_ModelAverage(dat,recurShock,respV,p,hmax,options)
% Auxiliary function for h-step impulse responses computed by a weighted average of different VAR submodels
% (modified based on Bruce Hansen's "cvar_ir.m" code, https://www.ssc.wisc.edu/~bhansen/progs/var.html)

% Inputs:
%	dat 	    nxm data matrix
%	p			max VAR lag order p >= 1
%	hmax		max horizon, hmax >= 1
%   recurShock  which shock?
%   respV       which respones?
%   options     numerical options for quadprog

% Outputs:
%	combination_irf		(hmax+1)    averaged impulse responses, horizon 0 to hmax
%	weights			    M x hmax 	model weights, by horizon (1 to hmax) and response variable
%   H                   m x m       choleskey decomp. of residual cov-var matrix


%% PREPARE

% Dimensions
n = size(dat,1)-p; % sample size
m = size(dat,2); % number of variables
k = m*p+1; % number of regressors

% Submodels: ordered as AR(1) to AR(p) and VAR(1) to VAR(p)
M = 2*p;
submodel_irf = zeros(hmax+1, M); % placeholder for IRF estimated by each submodel
combination_irf = zeros(hmax+1, 1); % placeholder for IRF averaged across submodels


%% FULL MODEL ESTIMATION

% Least-squares VAR
[~,By,sigma,Q,e,Bt,y,x] = VAR(dat,p);
H = chol(sigma,'lower'); % Warning: correspond to matrix C in our paper

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
jacobs_p = zeros(m*p,k*m);
G0 = zeros(k*m,hmax);
for h = 1:hmax
    lmax = min(h,p);
    the_irf_p = zeros(m,p);
    the_irf_p(:,1:lmax) = irf(:,h:-1:h-lmax+1);
    the_jacob = kron(eye(m), [0 the_irf_p(:)']) ...
                + Bt' * [zeros(1,k*m); jacobs_p];
    G0(:,h) = the_jacob(respV,:); % keep only response of interest
    jacobs_p = [the_jacob; jacobs_p(1:end-m,:)];
end


%% ESTIMATE SUBMODELS AND THEIR LOSS

% Auxiliary matrices in loss function
V_G0 = V*G0;
WI_G0 = reshape(Q\reshape(G0, [k m*hmax]), [k*m hmax]); % Peter J. Acklam, "MATLAB array manipulation tips and tricks", section 10.1.8

% Loop over submodels
Kr = zeros(M,hmax);

for r = 1:M-1
    
    % Identify constrained parameters and estimate constrained model
    B_sel = [true(m,1) false(m,m*p)]; % Default: only intercepts are unconstrained
    Btr = zeros(k,m);
    if r<=p % AR(1) to AR(p)
        B_sel(:,2:1+m*r) = repmat(eye(m)==1,1,r);
        for i=1:m
            [~,~,~,~,~,the_Bt_i] = VAR(dat(p+1-r:end,i),r);
            Btr([1 1+i+m*(0:r-1)],i) = the_Bt_i;
        end
    else % VAR(1) to VAR(p-1)
        r1 = r-p;
        B_sel(:,2:1+m*r1) = true;
        [~,~,~,~,~,Btr(1:1+r1*m,:)] = VAR(dat(p+1-r1:end,:),r1);
    end
    B_sel_t = B_sel';
    sel_vec = ~B_sel_t(:); % Constrained indices in theta=vec(B')
    
    % IRF
    er = y - x*Btr;
    sigmar = (er'*er) / (n-k+sum(sel_vec)/m);
    Hr = chol(sigmar,'lower');
    Byr = reshape(Btr(2:end,:),[m,p,m]);
    Byr = permute(Byr,[3,1,2]);
    irfr = IRF_SVAR(Byr,Hr(:,recurShock),hmax);
    submodel_irf(:,r) = irfr(respV,:);
    
    % Loss
    WI_sel_inv_V_G0 = WI(sel_vec,sel_vec)\V_G0(sel_vec,:);
    for h = 1:hmax
        Kr(r,h) = WI_G0(sel_vec,h)'*WI_sel_inv_V_G0(:,h);
    end
    
end


%% COMPUTE WEIGHTS FOR EACH HORIZONS

ir_diff = submodel_irf(2:end,:)'-var_ir;
weights = zeros(M, hmax);
ub = ones(M,1);
lb = zeros(M,1);
for h = 1:hmax
    ird = ir_diff(:,h);
    J0 = n*(ird*ird');
    K = - Kr(:,h)';
    w = quadprog(J0,K,[],[],ub',1,lb,ub,[],options); % compute optimal weights for submodels at each horizon h
    weights(:,h) = w;
    combination_irf(h+1,1) = submodel_irf(h+1,:)*w; % compute averaged IRF
end

end

