function varout = varest(y,var_par,smpl_par,levels)

% Computes VAR and covariance matrix of estimated parameters
%    
% Input:
%   y = txn matrix of series
%   var_par:
%     (1) nlag = number of VAR lags
%     (2) icomp = 0 Do not compute Companion Matrices
%                 1 Compute companion matrix ignoring contsant term
%                 2 Compute Companion matrix adding constant as last element of state -- loadings are zero if iconst = 0
%     (3) iconst  = 0 Do not include constant
%                   1 Include constant   		   
%     Note: icomp = 2 ... constant term is added to companion state vector (if iconst == 1)
%   smpl_par: sampling parameters
%      (1) calvec: Calendar vector
%      (2) nper:   number of periods per year
%      (3) nfirst: First observation for estimation
%      (4) nlast:  Last observation for estimation
%    		            
% Output:
%   varout
%      (1) betahat= estimated VAR coefficients (including constant if i_const == 1)
%                   each column has coefficients for a single equation
%      (2) seps = innovation covariance matrix
%      (3) resid =  VAR residuals
%      (4) coef: Companion form matrices
%           (a) Q
%           (b) M
%           (c) G
%           Companion Form Parameter, computed in icomp = 1;
%           Q, M, G for model written in SS form
%           y(t) = Q*z(t)
%           z(t) = M*z(t-1) + G*u(t)
%           var(u(t)) = I

% Extract parameters
nlag   = var_par.nlag;
icomp  = var_par.icomp;
iconst = var_par.iconst;
% Set Up VAR
ns = size(y,2);
T = size(y,1);
x = [ones(T,1),lagmatrix(y,1:nlag)];
if iconst == 0;
 x = x(:,2:end);  % Eliminate Constant if i_const == 0;
end;

[istart, iend] = smpl_HO(smpl_par);
istart = max(istart,1);
iend = min(iend,size(y,1));
trnd = (1:1:T)';
trnd = trnd(istart:iend);
y = y(istart:iend,:);
x = x(istart:iend,:);
resid = NaN*zeros(T,ns);

if levels % If data is in levels, estimate VECM

    % Johansen cointegration test
    [h,pValue,stat,cValue,mles] = jcitest(y,'Model','H*','Lags',nlag-1,'Test','maxeig');
    coint_rank = find(h{1,1:ns}==0,1)-1; % Estimated cointegration rank based on sequential tests
    mle = mles{1,coint_rank+1}; % MLE of the selected model
    varout.vecm.pval = pValue;
    varout.vecm.mle = mle;

    % Collect VAR parameters
    B_cell = cell(0);
    for l=1:nlag-1
        B_cell = [B_cell {mle.paramVals.(sprintf('B%d', l))}];
    end
    betahat_cell = vec2var(B_cell,mle.paramVals.A*mle.paramVals.B');
    betahat = [zeros(ns,1) cell2mat(betahat_cell)]'; % Add zeros for constant, since it will be ignored below
    seps = mle.EstCov;
    resid(istart+nlag:iend,:) = mle.res;

else % If data is differenced, do OLS
    
    tmp = packr([y x trnd]);
    y = tmp(:,1:ns);
    x = tmp(:,ns+1:end-1);
    trnd = tmp(:,end);
    betahat = x\y;
    e = y - x*betahat;
    ndf=size(x,1)-size(x,2);
    seps=(e'*e)/ndf;
    resid(trnd,:)=e;

end

if icomp > 0;
   %  
   %    Transform the VAR so that it is written in Standard form as:
   %    s(t)=P1*s(t-1) + P2*s(t-2) + ... + Pvarlag*s(t-varlag) + eps(t)
   %  

   % ---- Calculate Companion Matrix ---- ;
   if iconst == 0
      b = betahat';
      const_coef = zeros(ns,1);
   elseif iconst == 1
      b=betahat(2:end,:)';
      const_coef = betahat(1,:)';      % Coefficients on Constant Term ;
   else
      error('Invalid value of i_const in VAREST');
   end;
   comp=zeros(size(b,2),size(b,2));
   comp(1:size(b,1),:)=b;
   if size(b,2) > size(b,1);
     comp(size(b,1)+1:end,1:end-size(b,1))=eye(size(b,2)-size(b,1));
   end;
   %   -- Write Model in SS Form --
   %      y(t) = Q*z(t)
   %      z(t) = M*z(t-1) + G*u(t)
   %      var(u(t)) = I
   %  
   M=comp;
   Q=zeros(ns,size(M,1));
   Q(1:ns,1:ns)=eye(ns);
   G=zeros(size(M,1),ns);
   G(1:ns,1:ns)=(chol(seps))';
end;

if icomp == 2;  % Add Constant Term 
 G=[G ; zeros(1,size(G,2))];
 Q=[Q zeros(size(Q,1),1)];
 M=[M  zeros(size(M,1),1)];
 M=[M ; zeros(1,size(M,2))];
 M(end,end)=1;
 M(1:ns,end)=const_coef;
end;

% SAVE OUTPUT
varout.betahat = betahat;
varout.seps = seps;
varout.resid = resid;
varout.coef.Q = Q;
varout.coef.M = M;
varout.coef.G = G;



end