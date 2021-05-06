clear all;
addpath('Estimation_Routines');

% Test Hausman test code


%% Settings

% DGP (first variable is shock itself)
C = [1 0 0 0; 1 1 0 0; 1 0 1 0; 1 1 1 1]; % Impact impulse response matrix
mdl = varm('Constant', [0; ones(3,1)], ...
           'AR', {[0 0 0 0; 0 0.7 -0.2 0.1; 0 0 0.5 -0.3; 0 0 0 0.7], [0 0 0 0; 0 0.2 0.2 0; 0 0 0.3 0; 0 0 0 0]}, ...
           'Covariance', C*C'); % True VAR model

% Estimation settings
response_var = 3;       % Index of response variable
maxhorz = 8;            % Max impulse response horizon
p_estim = 1:2;          % Array of lag lengths used for estimation
signif_level = 0.1;     % Significance level

% Simulation settings
T = 1e4; % Sample size
numrep = 1e3; % No. of Monte Carlo repetitions
rng(20210506);

       
%% Model and true IRF

summarize(mdl);
irf_true_orthog = irf(mdl, 'NumObs', 1+maxhorz);
irf_true = irf_true_orthog(:,:,response_var)*(chol(C*C','lower')\C(:,1));


%% Simulate

num_p = length(p_estim);
hausman_stats = nan(numrep,num_p);
hausman_pvals = nan(numrep,num_p);
irf_tstats_VAR = nan(numrep,1+maxhorz,num_p);
irf_tstats_LP = nan(numrep,1+maxhorz,num_p);

timer = tic;

parfor i=1:numrep
    
    % Simulate data
    Y_sim = simulate(mdl, T);
    
    the_hausman_stats = nan(1,num_p);
    the_hausman_pvals = nan(1,num_p);
    the_irf_tstats_VAR = nan(1+maxhorz,num_p);
    the_irf_tstats_LP = nan(1+maxhorz,num_p);
    
    for j=1:num_p % For each estimation lag length...
        
        % Estimate VAR
        [~,the_By,the_Sigma,the_Sxx] = VAR(Y_sim,p_estim(j));
        the_C = chol(the_Sigma, 'lower');
        the_irf_VAR_all = IRF_SVAR(the_By,the_C(:,1),maxhorz);
        the_irf_VAR = the_irf_VAR_all(response_var,:)'/the_irf_VAR_all(1,1);
        
        % Hausman test
        [the_hausman_stats(j), the_hausman_pvals(j), the_irf_LP, the_varcov_VAR, the_varcov_LP] = ...
            IRF_Hausman(Y_sim,response_var,the_irf_VAR,the_By,the_Sigma,the_Sxx);
        
        % Test true IRF
        the_irf_VAR_diff = the_irf_VAR-irf_true;
        the_irf_tstats_VAR(:,j) = the_irf_VAR_diff./sqrt(diag(the_varcov_VAR));
        the_irf_LP_diff = the_irf_LP-irf_true;
        the_irf_tstats_LP(:,j) = the_irf_LP_diff./sqrt(diag(the_varcov_LP));
        
    end

    % Store results
    hausman_stats(i,:) = the_hausman_stats;
    hausman_pvals(i,:) = the_hausman_pvals;
    irf_tstats_VAR(i,:,:) = the_irf_tstats_VAR;
    irf_tstats_LP(i,:,:) = the_irf_tstats_LP;
    
    % Print progress
    if mod(i,ceil(numrep/50))==0
        fprintf('%3d%s\n', round(100*i/numrep), '%');
    end
    
end

elapsed_time = toc(timer);
disp('Elapsed min.:');
disp(elapsed_time/60);


%% Display results

% Expect rejection rate~=signif_level for p_estim < true lag length
% Expect rejection rate=signif_level for p_estim >= true lag length

cv = norminv(1-signif_level/2);

disp('Rejection rate: Hausman');
disp(mean(hausman_pvals<signif_level)); % p_estim along columns

disp('Rejection rate: VAR t-stat');
disp(squeeze(mean(abs(irf_tstats_VAR)>cv,1))); % Horizon along rows, p_estim along columns

disp('Rejection rate: LP t-stat');
disp(squeeze(mean(abs(irf_tstats_LP)>cv,1)));
