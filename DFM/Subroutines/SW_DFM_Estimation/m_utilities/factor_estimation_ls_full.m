function out = factor_estimation_ls_full(data, inclcode, est_par, levels)
% 2-28-2016, mww and po
% MODIFIFED 2/10/2018 to save residuals from UAR estimation 
%
% MODEL:  Y(t) = lambda*f(t) + u(t)
%         where Y(t) is ns x 1
%         f(t) is nfac x 1
%         This program uses the restricted data set to estimate the factors
%         and then uses the full dataset to estimate lambda's, AR model for
%         u and various summary statistics for the full dataset.
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
%     C. uniqueness
%        1. n_uarlag:  number or AR lags for uniquesess
% -- OUTPUT --
% out.est_data: data used to estimate facotors
% out.fac: output from factor_estimation_ls.m
% out.lam_mat: estimated lambda matrix (factor loadings)
% out.uar_coef_mat: univariate AR coefficients for uniquenesses;
% out.uar_ser_mat  = uar_ser_mat;
% out.uar_ser_mat: Estimated standard deviation for ;
% out.varout: VAR parameters for Factors -- see varest.m for description;
% out.r2: R-squared from regression of variables onto factors ;

% fac_est: Least squares estimates of factors
% lambda: Lease squares estimates of factor loadings
% tss: total sum of squares (standardized data)
% ssr: sum of squared residuals (standardized data)
% r2vec:  R2 for each series
% nobs: number of observations used for estimation
% nt: number of time series observations
% ns: number of cross-section observations


% PRELIMINARIES
n_series = size(data,2);                    % number of series
nfirst   = est_par.smpl_par.nfirst;         % start date
nlast    = est_par.smpl_par.nlast;          % end date
calvec   = est_par.smpl_par.calvec;         % calendar
nper     = est_par.smpl_par.nper;           % number of periods a year
n_uarlag = est_par.n_uarlag;                % number of AR lags
ntmin    = est_par.lambda.nt_min;           % minimum number of Obs

% USE SUBSET OF DATA TO ESTIMATE FACTORS
est_data = data(:,inclcode==1);
if levels
    est_data = [nan(1,size(est_data,2)); diff(est_data,1,1)]; % If data is in levels, estimate factors off first differences
end
lsout = factor_estimation_ls(est_data, est_par);
if levels
    lsout.fac_diff = lsout.fac;
    lsout.fac = cumsum_nan(lsout.fac); % If data was differenced, cumulate factors
end


% Compute estimates of factor loadings;
n_lc = 0;       % number of constraints placed on lambda
lambda_constraints_full = est_par.fac_par.lambda_constraints_full;
if size(lambda_constraints_full,2) > 1;
    lam_c_index = lambda_constraints_full(:,1);        % Which row of lambda: Constraints are then R*lambda = r
    lam_c_R     = lambda_constraints_full(:,2:end-1);  % R matrix
    lam_c_r     = lambda_constraints_full(:,end);      % r value
    n_lc = size(lambda_constraints_full,1);
end;
  
lam_mat = NaN(n_series,est_par.fac_par.nfac.total);
ismpl = smpl(calvec,nfirst,nlast,nper);
uar_coef_mat = NaN(n_series,n_uarlag);
uar_ser_mat = NaN(n_series,1);
uar_resid_mat = NaN(n_series,size(calvec,1));
r2_mat = NaN(n_series,1);                              % R-squared value
trend_tmp = (1:1:size(calvec,1))';
for is = 1:n_series;
    tmp = packr([data(ismpl==1,is) lsout.fac(ismpl==1,:) trend_tmp(ismpl==1)]);
    itmp = tmp(:,end);
    tmp = tmp(:,1:end-1+levels); % Include time trend if data is in levels
    if size(tmp,1) >= ntmin;
       y = tmp(:,1);
       x = [tmp(:,2:end), ones(size(tmp,1),1)];
       xxi = inv(x'*x);
       bols = xxi*(x'*y);
       b = bols;
       % Check for restrictions and impose;
       if n_lc > 0;
        ii = lam_c_index == is;
        if sum(ii) > 0;
            R = [lam_c_R(ii==1,:), zeros(sum(ii),1)];   % No constraints on constant term
            r = lam_c_r(ii==1,:);
            tmp1 = xxi*R';
            tmp2 = inv(R*tmp1);
            b = bols - tmp1*tmp2*(R*bols-r);
        end;
       end;
       lam_mat(is,:) = b(1:end-1-levels)';
       u = y - x*b;
       % Compute R-squared 
       ssr = sum(u.^2);
       ym = y - mean(y);
       tss = sum(ym.^2);
       r2_mat(is) = 1-(ssr/tss);
       % Compute AR model for errors
       if r2_mat(is) < 0.9999;
        [arcoef, ser, ar_resid] = uar(u,n_uarlag);    % AR Coefficients and ser 
       else;
        arcoef = zeros(n_uarlag,1);
        ser = 0.0;
        ar_resid = NaN(size(u,1),1);
       end;
       uar_coef_mat(is,:) = arcoef'; 
       uar_ser_mat(is,1) = ser;
       uar_resid_mat(is,itmp) = ar_resid';
    end;
end;

varout = varest(lsout.fac,est_par.var_par,est_par.smpl_par,levels);

% SAVE OUTPUT
out.est_data      = est_data;
out.fac           = lsout.fac;
out.lam_mat       = lam_mat;
out.uar_coef_mat  = uar_coef_mat;
out.uar_ser_mat   = uar_ser_mat;
out.uar_resid_mat = uar_resid_mat;
out.varout        = varout;
out.r2            = r2_mat;

if levels
    out.vecm = varout.vecm;
end

end