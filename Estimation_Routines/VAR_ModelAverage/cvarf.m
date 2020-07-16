function [combination_forecast,weights,combination_beta] = cvarf(dat,p,hmax)

% Creates h-step point forecasts from VAR(p) model with intercept

% Inputs:
%	dat 			nxm data matrix
%	p			VAR order p >= 1
%	hmax			forecast horizon, hmax >= 1

% Outputs:
%	combination_forecast	hmax x m	point forecasts, for horizons h=1,...,hmax, for each variable
%	weights			M x hmax x m	model weights, for forecast variable and forecast horizon (M=2p)
%	combination_beta	k x hmax x m	forecast coefficient vectors, by forecast variable and horizon (k=mp+1)

n = size(dat,1)-p;
m = size(dat,2);
k = m*p+1;
y = dat(1+p:n+p,:);
x = ones(n,k);
xf = ones(1,k);
for i = 1:p
    x(:,(i-1)*m+1:i*m) = dat(1+p-i:n+p-i,:);
    xf(1,(i-1)*m+1:i*m) = dat(1+n+p-i,:);
end

% Full Model Estimation
Bt = (x'*x)\(x'*y);
B = Bt';
theta = Bt(:);
Q = (x'*x) / n;
C = inv(chol(Q));
Qinv = C*(C');
e = y - x*B';
xe = zeros(n,k*m);
for i = 1:m
  xe(:,(i-1)*k+1:i*k) = x .* (e(:,i)*ones(1,k));
end
omega = (xe'*xe) / (n-k);
WI = kron(eye(m),Qinv);
V = WI*omega*WI;
J = [eye(m*(p-1)), zeros(m*(p-1),m+1); zeros(1,m*p), 1];
P = [B; J];

G = zeros(k*m,k,m,hmax);
for j = 1:m
  for h = 1:hmax
    Gjh = zeros(k*m,k);
    for i = 1:h
      tmp = P^(h-i);
      Gjh = Gjh + kron(tmp(j,1:m)',P^(i-1));
    end
    G(:,:,j,h) = Gjh;
  end
end

M = 2*p;
submodel_betas = zeros(k, M, hmax, m); 
submodel_forecasts = zeros(hmax, M, m);
Kr = zeros(M,m,hmax);

for r = 1:M		% Submodels
  if r<=p
    R = zeros(0,0);
    R2 = [zeros(m*r,m*(p-r));eye(m*(p-r));zeros(1,m*(p-r))];
    for i = 1:m
      Ci = eye(m);
      Ci(:,i) = [];
      R1 = [kron(eye(r),Ci); zeros(m*(p-r)+1,(m-1)*r)];
      R = blkdiag(R,[R1,R2]);
    end
    WR = (WI*R)*((R'*WI*R)\(R'));
  elseif r<(2*p)
    r1 = r-p;
    R2 = [zeros(m*r1,m*(p-r1));eye(m*(p-r1));zeros(1,m*(p-r1))];
    R = kron(eye(m),R2);
    WR = (WI*R)*((R'*WI*R)\(R'));
  elseif r==(2*p)
    WR = zeros(k*m);
  end
  thetar = theta - WR*theta;
  Br = (reshape(thetar,k,m))';
  Pr = [Br; J];
  Ph = eye(k);
  for h = 1:hmax
    Ph = Ph*Pr;
    betar = Ph(1:m,:)';			% Forecast Coefficients
    submodel_betas(:,r,h,:) = betar;
    submodel_forecasts(h,r,:) = (xf*betar)';	% Point Forecasts
    for j = 1:m
      Gj  = G(:,:,j,h);
      Kr(r,j,h) = trace(Q*Gj'*WR*V*Gj);
    end
  end
end

combination_beta = zeros(k, hmax, m); 
combination_forecast = zeros(hmax, m);
weights = zeros(M, hmax, m);

ub = ones(M,1);
lb = zeros(M,1);
options = optimoptions('quadprog','Display','off');
for i = 1:m
  for h = 1:hmax
    betadiff = submodel_betas(:,:,h,i) - submodel_betas(:,end,h,i)*ones(1,M);
    J = n*betadiff'*Q*betadiff;
    J = (J+J')/2;
    K = - Kr(:,i,h)';
    w = quadprog(J,K,[],[],ub',1,lb,ub,[],options);
    weights(:,h,i) = w;
    combination_beta(:,h,i) = submodel_betas(:,:,h,i)*w;
    combination_forecast(h,i) = submodel_forecasts(h,:,i)*w;
  end
end


