function DFM = DFM_est(n_factors,n_lags_fac,n_lags_uar, reorder, levels);
% Function for estimating parameters in the encompassing DFM model
% (Revised based on Mark Watson's MATLAB script)

%% PREPARATIONS

small = 1.0e-10;
big = 1.0e+6;

% ----------- Sample Period, Calendars and so forth
[dnobs_m,calvec_m,calds_m] = calendar_make([1959 1],[2014 12],12);  % Monthly Calendar
[dnobs_q,calvec_q,calds_q] = calendar_make([1959 1],[2014 4],4);    % Quarterly Calendar

% -- Load Data
load_data=1;
  % Demeaning Parameters
  i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;

  bw_bw = 100;   % Bi-Weight Parameter for local demeaning
datain_all;      % datain_all reads in the full dataset .. all variables, etc. saved in datain.xx

%% ESTIMATION

% Factor Parameters
n_fac = n_factors;
est_par.fac_par.nfac.unobserved = n_fac;
est_par.fac_par.nfac.observed = 0;
est_par.fac_par.nfac.total = n_fac;

% Sampling parameters
est_par.smpl_par.nfirst = [1959 3];       % start date
est_par.smpl_par.nlast  = [2014 4];       % end date
est_par.smpl_par.calvec = datain.calvec;  % calendar
est_par.smpl_par.nper   = 4;              % number of periods a year

% Factor analysis parameters
est_par.fac_par.nt_min                  = 20;     % min number of obs for any series used to est factors
est_par.lambda.nt_min                   = 40;     % min number of obs for any series used to estimate lamba, irfs, etc.
est_par.fac_par.tol                     = 10^-8;  % precision of factor estimation (scaled by by n*t)

% Restrictions on factor loadings to identify factors
est_par.fac_par.lambda_constraints_est  = 1;  % no constraints on lambda
est_par.fac_par.lambda_constraints_full = 1;  % no constraints on lambda

% VAR parameters for factors
est_par.var_par.nlag   = n_lags_fac;    % number of lags
est_par.var_par.iconst = 1;    % include constant
est_par.var_par.icomp  = 1;    % compute companion form of model .. excluding constant

% yit equation parameters
est_par.n_uarlag = n_lags_uar;  % number of arlags for uniqueness

% Matrices for storing results
n_series = size(datain.bpdata,2);

fac_est_out = factor_estimation_ls_full(datain.bpdata,datain.bpinclcode,est_par,levels);                  % estimation

% Save some output;
% Calendar, etc;
factor_model.calvec = est_par.smpl_par.calvec;
% Factors
factor_model.nfac = n_fac; % number of factors
factor_model.fac = fac_est_out.fac;   % Factor Estimates
% Variable indentifiers
factor_model.bpnamevec = datain.bpnamevec;
factor_model.bplabvec_long = datain.bplabvec_long;
factor_model.bplabvec_short = datain.bplabvec_short;
factor_model.bptcodevec = datain.bptcodevec;
% Factor VAR
factor_model.var_nlag = est_par.var_par.nlag;
factor_model.var_M = fac_est_out.varout.coef.M;
factor_model.var_Q = fac_est_out.varout.coef.Q;
factor_model.var_G = fac_est_out.varout.coef.G;
factor_model.var_seps = fac_est_out.varout.seps;
factor_model.var_resid = fac_est_out.varout.resid;
% Univariate AR
factor_model.uar_coef_mat = fac_est_out.uar_coef_mat;
factor_model.uar_ser_mat = fac_est_out.uar_ser_mat;
factor_model.uar_resid_mat = fac_est_out.uar_resid_mat;
% Factor loadings
factor_model.lam_mat = fac_est_out.lam_mat;

%% COLLECT FINAL MATRICES

DFM.Lambda     = fac_est_out.lam_mat;
DFM.Phi        = fac_est_out.varout.coef.M;
DFM.n_lags_fac = size(DFM.Phi,1) / n_factors;
DFM.Sigma_eta  = fac_est_out.varout.seps; % var-cov matrix for reduced-form shocks
DFM.sigma_v    = factor_model.uar_ser_mat;
DFM.delta      = factor_model.uar_coef_mat;
DFM.v          = fac_est_out.res_mat;

if levels
    DFM.vecm = fac_est_out.vecm;
end

DFM.bpnamevec = datain.bpnamevec; % variable name
DFM.bplabvec_long = datain.bplabvec_long;
DFM.bplabvec_short = datain.bplabvec_short;
DFM.bptcodevec = datain.bptcodevec; % transformation code

DFM.fac = fac_est_out.fac;
DFM.fac_res = fac_est_out.varout.resid;
DFM.r2 = fac_est_out.r2;
DFM.factor_shock_time_range = [1959 2014.75 4]; % (start, end, freq)

%% REORDER VARIABLES TO MATCH VARIABLE LIST IN S&W (2016)

DFM.Lambda(reorder, :) = DFM.Lambda;
DFM.sigma_v(reorder) = DFM.sigma_v;
DFM.delta(reorder, :) = DFM.delta;
DFM.v(reorder, :) = DFM.v;

DFM.bpnamevec(reorder) = DFM.bpnamevec;
DFM.bplabvec_long(reorder) = DFM.bplabvec_long;
DFM.bplabvec_short(reorder) = DFM.bplabvec_short;
DFM.bptcodevec(reorder) = DFM.bptcodevec;

DFM.r2(reorder) = DFM.r2;

end