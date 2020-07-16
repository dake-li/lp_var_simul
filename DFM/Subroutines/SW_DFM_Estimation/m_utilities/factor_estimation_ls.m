function lsout = factor_estimation_ls(est_data, est_par)
% 2-26-2016, mww and po
%
% MODEL:  Y(t) = lambda*f(t) + u(t)
%         where Y(t) is ns x 1
%         f(t) is nfac x 1
%
% -- INPUT --
% est_data: Data used for estimating factors  (nt x ns)
% est_par: estimation parameters
%     A. smpl_par: sampling parameters
%         1. calvec: Calendar vector
%         2. nper:   number of periods per year
%         3. nfirst: First observation for estimation
%         4. nlast:  Last observation for estimation
%     B. fac_par: factor estimation parameters
%         1. w: A set of observed factors. Not used if no observed factors
%         2. lambda_constraints: These are linear constraints on the lambda parameters
%                                Each row of lambda_constraints corresponds to a single constraint
%                                .. first column: "i" ... which row of lambda is constrained;
%                                .. call this row lam(i)
%                                The constraint is then R*lam(i)' = r
%                                R is given in cols 2-nfac+1 of lambda_constraints
%                                r is given in last column of lambda_constraints
%                                set = 1 if not used;
%         3. ntmin:  Minimum number of observations for any series used to estimate factors
%         4. tol: Precision of estimate
%                 Loop terminates when objective function changes by less than tol/(nt*ns)
%         5. nfac:
%            (a) unobserved: number of unobserved factors
%            (b) observed: number of observed factors
%            (c) total: total number of factors 
%   
% -- OUTPUT --
% fac_est: Least squares estimates of factors
% lambda: Lease squares estimates of factor loadings
% tss: total sum of squares (standardized data)
% ssr: sum of squared residuals (standardized data)
% r2vec:  R2 for each series
% nobs: number of observations used for estimation
% nt: number of time series observations
% ns: number of cross-section observations


% Preliminaries
% extract estimation parameters
smpl_par           = est_par.smpl_par;
lambda_constraints = est_par.fac_par.lambda_constraints_est;
nt_min             = est_par.fac_par.nt_min;
tol                = est_par.fac_par.tol;
nfac_u             = est_par.fac_par.nfac.unobserved;
nfac_o             = est_par.fac_par.nfac.observed;
nfac_t             = est_par.fac_par.nfac.total;
if nfac_o > 0;
  w = est_par.fac_par.w;
end;

% Estimate factors with unbalanced panel by LS -- standardize data first
  % Sample period
  [istart, iend] = smpl_HO(smpl_par);
  istart = max(istart,1);
  iend = min(iend,size(est_data,1));
  
  % Estimate Factors 
  xdata = est_data(istart:iend,:);
  nt = size(xdata,1);
  ns = size(xdata,2);
   
  % Mean and Standard Deviation
  xmean = nanmean(xdata)';                                  % mean (ignoring NaN)
  mult = sqrt((sum(~isnan(xdata))-1)./sum(~isnan(xdata)));  % num of non-NaN entries for each series
  xstd = (nanstd(xdata).*mult)';                            % std (ignoring NaN)
  xdata_std = (xdata - repmat(xmean',nt,1))./repmat(xstd',nt,1);   % standardized data
  
  if nfac_o > 0;
      wdata = w(istart:iend,:);
      if sum(sum(isnan(wdata))) ~= 0;
          error('w contains missing values over sample period, processing stops');
      end;
  end;
  
  n_lc = 0;            % number of constrains placed on lambda
  if size(lambda_constraints,2) > 1;
     lam_c_index = lambda_constraints(:,1);     % Which row of lambda: Constraints are then R*lambda = r
     lam_c_R = lambda_constraints(:,2:end-1);   % R matrix
     lam_c_r = lambda_constraints(:,end);       % r value
     lam_c_r_scl = lam_c_r./xstd(lam_c_index);  % Adjusted for scaling
     n_lc = size(lambda_constraints,1);
  end;
     
  % Compute Total Sum of Squares
  tss = 0;
  nobs = 0;
  for is = 1:ns;
      tmp = xdata_std(:,is);     % select series
      tmp = tmp(isnan(tmp)==0);  % drop NaN
      tss = tss+sum(tmp.^2);     % add to tss
      nobs = nobs+size(tmp,1);   % add to n*T
  end;
  
  % Estimate factors using balanced panel
  if nfac_u > 0;
   xbal = packr(xdata_std')';
   %[coef,score,latent]=princomp(xbal);
   [coef,score,latent] = pca(xbal);
   f = score(:,1:nfac_u);
   fa = f;
   if nfac_o > 0;
      fa = [wdata f];
   end;
   lambda = NaN*zeros(ns,nfac_t);
  else;
   fa = wdata;
  end;
  
  diff = 100;
  ssr = 0;
  while diff>tol*(nt*ns)
      ssr_old = ssr;
    for i = 1:ns; 
        tmp=packr([xdata_std(:,i) fa]);
        if size(tmp,1) >= nt_min;
  	       y=tmp(:,1);
     	     x=tmp(:,2:end);
           xxi = inv(x'*x);
           bols = xxi*(x'*y);
           b = bols;
           if n_lc > 0;
             % Check for restrictions and impose;
             ii = lam_c_index == i;
             if sum(ii) > 0;
                 R = lam_c_R(ii==1,:);
                 r_scl = lam_c_r_scl(ii==1,:);
                 tmp1 = xxi*R';
                 tmp2 = inv(R*tmp1);
                 b = bols - tmp1*tmp2*(R*bols-r_scl);
             end;
           end;
           lambda(i,:)= b';
        end;
    end;
    edata = xdata_std;
    if nfac_u > 0;
       if nfac_o > 0;
          edata = xdata_std - fa(:,1:nfac_o)*lambda(:,1:nfac_o)';
       end;
       for t = 1:nt;
          tmp=packr([edata(t,:)' lambda(:,nfac_o+1:end)]);
          y=tmp(:,1);
     	    x=tmp(:,2:end);
          b = x\y;
          f(t,:) = b';
       end; 
       fa = f;
       if nfac_o > 0;
         fa = [wdata f];
       end;
    end;
    % Compute residuals 
    e = xdata_std-fa*lambda';
    ssr = sum(nansum(e.^2));
    diff = abs(ssr_old - ssr);
  end;
  
  f_est = fa;
  
  % Compute R2 for each series
  r2vec = NaN*zeros(ns,1);
  for i = 1:ns;
      tmp=packr([xdata_std(:,i) f_est]);
      if size(tmp,1) >= nt_min;
  	    y=tmp(:,1);
        x=tmp(:,2:end);
        b = x\y;
        e = y -x*b;
        r2_ssr = e'*e;
        r2_tss = y'*y;
        r2vec(i) = 1-r2_ssr/r2_tss;
      end;
  end;
  
  fac_est = NaN(size(est_data,1),nfac_t);
  fac_est(istart:iend,:)=f_est;
  lambda = lambda.*repmat(xstd,1,nfac_t);

  lsout.fac = fac_est;   lsout.lambda = lambda;
  lsout.tss = tss;       lsout.ssr = ssr;
  lsout.r2vec = r2vec;   lsout.nobs = nobs;
  lsout.nt = nt;         lsout.ns = ns;
  
end

