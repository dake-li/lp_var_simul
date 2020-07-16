% HAR factor model compute ACVs
%

clear all;
small = 1.0e-10;
big = 1.0e+6;

% -- File Directories  
figdir = 'fig/';
outdir = 'out/';
matdir = 'mat/';
matcvdir = 'mat_cv/';
procs_dir = '../m_utilities/';            % Directory for necessary matlab_procs and functions
p = path;   % Update path to access procs_dir                              
p1=path(procs_dir,p);

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);  % Monthly Calendar
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);    % Quarterly Calendar

% -- Load Data
load_data=0;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;
    
  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
% datain_all;      % datain_all reads in the full dataset .. all variables, etc. saved in datain.xx
fstr = [matdir 'har_fac_' num2str(i_demean)];load(fstr,'har_fac'); 

% Calendar, etc;
calvec = har_fac.calvec; % number of factors
fac = har_fac.fac;   % Factor Estimates
% Variable indentifiers
bpnamevec = har_fac.bpnamevec;
bplabvec_long = har_fac.bplabvec_long;
bplabvec_short = har_fac.bplabvec_short;
bptcodevec = har_fac.bptcodevec;
% Factor VAR
var_nlag = har_fac.var_nlag;
M = har_fac.var_M;
Q = har_fac.var_Q;
G = har_fac.var_G;
%           Q, M, G for model written in SS form
%           y(t) = Q*z(t)
%           z(t) = M*z(t-1) + G*u(t)
%           var(u(t)) = I


var_seps = har_fac.var_seps;
var_resid = har_fac.var_resid;
% Univariate AR
uar_coef_mat = har_fac.uar_coef_mat;
uar_ser_mat = har_fac.uar_ser_mat;
uar_resid_mat = har_fac.uar_resid_mat;
% Factor loadings
lam_mat = har_fac.lam_mat;

n_acv = 500;   % Number of autocovariances to compute
n_factor = size(Q,1);
n_series = size(lam_mat,1);

% Compute ACV of 
GG = G*G';
n_size = size(M,1);
tmp = eye(n_size^2) - kron(M,M);
tmp = inv(tmp);
vec_GG = reshape(GG,n_size^2,1);
vec_sig = tmp*vec_GG;
sig = reshape(vec_sig,n_size,n_size);
ACV_Factor = NaN(n_factor,n_factor,n_acv+1);
ACV_Factor(:,:,1) = Q*sig*Q';
for i = 1:n_acv;
    sig = M*sig;
    ACV_Factor(:,:,i+1) = Q*sig*Q';
end;

acv_mat = NaN(n_acv+1,n_series);
for jj = 1: n_series;
  lam = lam_mat(jj,:)';
  acv_f = NaN(n_acv+1,1);
  acv_f(1) = lam'*squeeze(ACV_Factor(:,:,1))*lam;
  for i = 1:n_acv;
    acv_f(1+i) = lam'*squeeze(ACV_Factor(:,:,1+i))*lam;
  end;
  ar_coef = uar_coef_mat(jj,:);
  u_var = uar_ser_mat(jj)^2;
  acv_u = arma_acv(ar_coef,0,n_acv );
  acv_u = u_var*acv_u;
  acv_t = acv_u+acv_f;
  acv_mat(:,jj) = acv_t;
end;

fstr = [matdir 'har_fac_acvmatrix'];save(fstr,'acv_mat'); 

path(p);  % Reset path