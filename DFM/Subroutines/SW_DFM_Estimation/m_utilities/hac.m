function [betahat,vbeta,se_beta,ser,rbarsq] = hac(y,x,nma,ikern);
%{
     Modified by MWW, 12-28-96

     Procedure for estimating the regression y = xbeta+ u
     The procedure produces the OLS estimate of b
     and a hetero/autocorrelation consistent estimate of
     the autocorrelation matrix.

Input:
     y = tx1
     x = txk
     nma=truncation parameter (nma=0, White SEs) 
     ikern = kernel indicator
             1 => triangular
             2 => rectangular

Output:
     Beta = OLS estimate of beta (kx1)
    VBeta = Robust estimate of covariance matrix of beta (kxk)
            (Note this is computed used PINV, if X'X is singular

%}
xx=x'*x;
xxi = inv(xx);
betahat=x\y;
u=y-x*betahat;

z = x.*repmat(u,1,size(x,2));
v=zeros(size(x,2),size(x,2));


% Form Kernel 
kern=zeros(nma+1,1);
for ii = 0:nma;
 kern(ii+1,1)=1;
 if nma > 0;
  if ikern == 1; 
    kern(ii+1,1)=(1-(ii/(nma+1))); 
  end;
 end;
end;

% Form Hetero-Serial Correlation Robust Covariance Matrix 
for ii = -nma:nma;
 if ii <= 0; 
    r1=1; 
    r2=size(z,1)+ii; 
 else; 
    r1=1+ii; 
    r2=size(z,1); 
 end;
 v=v + kern(abs(ii)+1,1)*(z(r1:r2,:)'*z(r1-ii:r2-ii,:));
end;

vbeta=xxi*v*xxi';

% Compute Other Statistics of interest
y_m = y - mean(y);
tss = y_m'*y_m;
ess = u'*u;
ndf = size(x,1)-size(x,2);
ser = sqrt(ess/ndf);
rbarsq = 1 - (ess/ndf)/(tss/(size(x,1)-1));
se_beta = sqrt(diag(vbeta));


end