function [ lm, lmr, lsbreak ] = qlra(y,x1,x2,ccut,nma);
% -- Computes QLR Statistic -- Robust and Non-Robust
%
%  Input:
%   y = lhv data
%   x1 = rhv data with fixed coefficients under null and alternative
%        (input a scalar [1,1] matrix if all variables
%         are allowed to vary)
%   x2 = rhv data with fixed coefficients under null and
%        time varying coefficients under alternative
%   ccut=endpoints for sequential chow regressions
%
%    nma -- number of MA components for HAC Matrix
%           (0 => White Hetero Robust)
%
%    Output:
%
%    LM -- QLR Statistic
%    LMR -- QLR Robust Statistic
%    lsbreak -- least squares estimate of break date 

 nobs=length(y); ktrim=floor(ccut*nobs); n1t=ktrim; n2t=nobs-ktrim;
 k=size(x2,2);

 big=1.0e+15;
 lr=-big*ones(nobs,1);
 lrr=-big*ones(nobs,1);
 ssrvec=big*ones(nobs,1);

 % N-W Kernel 
 kern=zeros(nma+1,1);
 for ii = 0:nma;
  kern(ii+1,1)=1;
  if nma > 0;
   kern(ii+1,1)=(1-(ii/(nma+1))); 
  end;
 end;
  
 x=x2; if length(x1) ~= 1; x=[x1,x2]; end;
% full sample @
  sxy=x'*y; sxx=x'*x; syy=y'*y;

  sx2y=x2(1:n1t-1,:)'*y(1:n1t-1,:);
  sx2x2=x2(1:n1t-1,:)'*x2(1:n1t-1,:);
  sxx2=x(1:n1t-1,:)'*x2(1:n1t-1,:);

  for i = n1t:n2t;
    sx2y=sx2y+x2(i,:)'*y(i,:);
    sx2x2=sx2x2+x2(i,:)'*x2(i,:);
    sxx2=sxx2+x(i,:)'*x2(i,:);
    mxx=[[sxx,sxx2];[sxx2',sx2x2]];
    mxy=[sxy;sx2y];
    mxxi=inv(mxx);
    bhat=mxxi*mxy;
    w=[x,[x2(1:i,:);zeros(nobs-i,size(x2,2))]];
    e1=y-w*bhat;
    ssr=e1'*e1;
    ssrvec(i)=ssr;
    s2=ssr/(length(e1)-length(bhat));

  % Form Hetero-Serial Correlation Robust Covariance Matrix @
    we=w.*repmat(e1,1,size(w,2));
    v=zeros(size(we,2),size(we,2));
    for ii =-nma:nma;
     if ii <= 0; 
      r1=1; 
      r2=length(we)+ii; 
     else; 
      r1=1+ii; 
      r2=length(we); 
     end;
     v=v + kern(abs(ii)+1,1)*(we(r1:r2,:)'*we(r1-ii:r2-ii,:));
    end;
 
    vbetar=mxxi*(v)*mxxi;
    vbeta=s2*mxxi;
    bet=bhat(size(x,2)+1:length(bhat),1);
    vbetar=vbetar(size(x,2)+1:length(bhat),size(x,2)+1:length(bhat));
    vbeta=vbeta(size(x,2)+1:length(bhat),size(x,2)+1:length(bhat));
    lr(i)=bet'*(inv(vbeta))*bet;
    lrr(i)=bet'*(inv(vbetar))*bet;
  end;
    
    lm=max(lr);
    lmr=max(lrr);
    [tmp,lsbreak]=min(ssrvec);

