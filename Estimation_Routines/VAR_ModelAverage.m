function [combination_irf,weights,H] = VAR_ModelAverage(dat,recurShock,respV,p,hmax)
% Creates h-step impulse response functions from estimated VAR(p) model with intercept
% (modified based on Bruce Hansen's "cvar_ir.m" code)
% https://www.ssc.wisc.edu/~bhansen/progs/var.html

% Inputs:
%	dat 			nxm data matrix
%	p			VAR order p >= 1
%	hmax			horizon, hmax >= 1
%   recurShock  which shock?
%   respV       which respones?

% Outputs:
%	combination_irf		(hmax+1) orthogonalized impulse responses, horizon 0 to hmax
%	weights			M x hmax 	model weights, by horizon (1 to hmax) and response variable
%   H               m x m       choleskey decomp. of residual cov-var matrix

% prepare
n = size(dat,1)-p;
m = size(dat,2);
k = m*p+1;
y = dat(1+p:n+p,:);
x = ones(n,k);
for i = 1:p
    x(:,(i-1)*m+1:i*m) = dat(1+p-i:n+p-i,:);
end

% Full Model Estimation
Bt = (x'*x)\(x'*y);
B = Bt';
theta = Bt(:);
Q = (x'*x) / n;
C = inv(chol(Q));
Qinv = C*(C');
e = y - x*B';
sigma = (e'*e) / (n-k);
H = chol(sigma,'lower');
xe = zeros(n,k*m);
for i = 1:m
  xe(:,(i-1)*k+1:i*k) = x .* (e(:,i)*ones(1,k));
end
omega = (xe'*xe) / (n-k);
WI = kron(eye(m),Qinv);
V = WI*omega*WI;
J = [eye(m*(p-1)), zeros(m*(p-1),m+1); zeros(1,m*p), 1];
P = [B(:,1:m*p), zeros(m,1); J];

M = 2*p;
submodel_irf = zeros(hmax+1, M); 
var_ir = zeros(1,hmax);
combination_irf = zeros(hmax+1, 1); 
weights = zeros(M, hmax);
Kr = zeros(M,hmax);
ir_diff = zeros(M,hmax);

submodel_irf(1,M) = H(respV, recurShock);
combination_irf(1,1) = H(respV, recurShock);

Ph0 = [H; zeros(k-m,m)];
Ph = Ph0;
for h = 1:hmax
  Ph = P*Ph;
  var_ir(1,h) = Ph(respV, recurShock);
  submodel_irf(h+1,M) = Ph(respV, recurShock);
end

% derivative and Avar
G0 = zeros(k*m,hmax);
VG = zeros(1,hmax);
for h = 1:hmax
    Gh = zeros(k*m,m^2);
    T2 = Ph0;
    for i = 1:h
      T1 = (P^(h-i))';
      Gh = Gh + kron(T1(1:m,1:m),T2);
      T2 = P*T2;
    end
    G0(:,h) = Gh(:, (respV - 1) * m + recurShock); % keep only response of interest
    % VG(:,:,h) = inv(Gh'*V*Gh); % asymp var as weighting matrix
    % VG(:,:,h) = eye(m^2); % identity weighting matrix
    VG(1,h) = 1; % trivial weighting matrix
end

for r = 1:M		% Submodels
  if r<=p % AR(1) to AR(p)
    R = zeros(0,0);
    R2 = [zeros(m*r,m*(p-r));eye(m*(p-r));zeros(1,m*(p-r))];
    for i = 1:m
      Ci = eye(m);
      Ci(:,i) = [];
      R1 = [kron(eye(r),Ci); zeros(m*(p-r)+1,(m-1)*r)];
      R = blkdiag(R,[R1,R2]);
    end
    WR = (WI*R)*((R'*WI*R)\(R'));
  elseif r<(2*p) % VAR(1) to VAR(p)
    r1 = r-p;
    R2 = [zeros(m*r1,m*(p-r1));eye(m*(p-r1));zeros(1,m*(p-r1))];
    R = kron(eye(m),R2);
    WR = (WI*R)*((R'*WI*R)\(R'));
  elseif r==(2*p)
    WR = zeros(k*m);
  end
  thetar = theta - WR*theta;
  Br = (reshape(thetar,k,m))';
  er = y - x*Br';
  sigmar = (er'*er) / (n-k+size(R,2)/m);
  Hr = chol(sigmar,'lower');
  Pr = [Br(:,1:(m*p)), zeros(m,1); J];
  Ph = [Hr; zeros(k-m,m)];
  submodel_irf(1,r) = Hr(respV, recurShock);
  for h = 1:hmax
    Ph = Pr*Ph;
    Phm = Ph(respV, recurShock);
    diff = Phm - var_ir(1,h);
    ir_diff(r,h) = diff;
    submodel_irf(h+1,r) = Phm;
    Gh = G0(:,h);
    Vg = VG(1,h);
    Kr(r,h) = trace(Vg*Gh'*WR*V*Gh);
  end
end

ub = ones(M,1);
lb = zeros(M,1);
options = optimoptions('quadprog','Display','off');
for h = 1:hmax
    ird = ir_diff(:,h);
    Vg = VG(1,h);
    J0 = n*ird*Vg*ird';
    J0 = (J0+J0')/2;
    K = - Kr(:,h)';
    w = quadprog(J0,K,[],[],ub',1,lb,ub,[],options);
    weights(:,h) = w;
    combination_irf(h+1,1) = submodel_irf(h+1,:)*w;
end

end

