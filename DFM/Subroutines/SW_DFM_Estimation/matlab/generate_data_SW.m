function [y_data,shocks_data] = generate_data_SW(DFM,settings)
% Function for generating simulated data using calibrated DFM

%%%%%%directly use Mark Watson's MATLAB script%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Generate variables from factor model using normal innovations
% clear all;
small = 1.0e-10;
big = 1.0e+6;
% rng(63876498);
% Demeaning Parameters .. model for data when estimates
i_demean = 1;  % 0 Do Nothing
               % 1 Eliminate low-frequency by local Demeaning
               % 2 Eliminate low-frequency trend by full-sample demeaning;


% -- File Directories  
% figdir = 'fig/';
% outdir = 'out/';
% matdir = 'mat/';
% sdata_dir = 'sdata/';
% procs_dir = '../m_utilities/';            % Directory for necessary matlab_procs and functions
% p = path;   % Update path to access procs_dir                              
% p1=path(procs_dir,p);
% fstr = [matdir 'har_fac_' num2str(i_demean)];load(fstr);

nT = settings.simul.T;
nT_burn = settings.simul.T_burn;

lam_mat = DFM.Lambda;
uar_coef_mat = DFM.delta;
uar_ser_mat = DFM.sigma_v;
M = DFM.Phi;

if DFM.n_lags_fac == 1
    Q = eye(settings.DFM.n_fac);
    G = chol(DFM.Sigma_eta)';
else
    Q = [eye(settings.DFM.n_fac),zeros(settings.DFM.n_fac,(DFM.n_lags_fac-1)*settings.DFM.n_fac)];
    G = [chol(DFM.Sigma_eta)'; zeros((DFM.n_lags_fac-1)*settings.DFM.n_fac,settings.DFM.n_fac)];
end

% Generate draws of factors
T = nT;           % Number of observations to generate
T_initial = nT_burn;   % Initial obs
nrep = 1;
% Generate Factors
[F_mat,shocks_data] = generate_data_VAR_companion(T,nrep,Q,M,G,T_initial);

ny = size(lam_mat,1);
y_data = NaN(T,nrep,ny);

for i = 1:ny;
    % i
    s_mat = NaN(T,nrep);
    lam = lam_mat(i,:);
    arcoef = uar_coef_mat(i,:);
    se_ar = uar_ser_mat(i);
    u_mat = generate_data_univariate_ar(T,nrep,arcoef,se_ar,T_initial);
    for t = 1:T;
        F = F_mat(t,:,:);
        F = permute(F,[2,3,1]);
        s_mat(t,:) = lam*F+u_mat(t,:);
    end;
    y_data(:,:,i) = s_mat;
    % fstr = [sdata_dir 'sdata_normal_' num2str(i_demean) '_' num2str(i)];save(fstr,'s_mat','-v7.3');
end;
y_data = permute(y_data,[1,3,2]); % T * v * rep


% path(p);  % Reset path

end

