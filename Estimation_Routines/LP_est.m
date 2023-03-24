function [IRF,n_lags_est] = LP_est(data_sim,settings,bias_corrected);
% Function for estimating IRFs using least-squares LP

% preparations

run('Estimation_Setup'); % common setup for all estimation methods

% estimate IRF via LP

[IRF_resp,w] = IRF_LP(Y,recursiveShock,responseV,nlags,IRF_hor - 1); % IRF to one unit of shock
if bias_corrected == 1
    IRF_resp = LP_CorrectBias(IRF_resp, w); % Herbst-Johanssen bias correction
end

% normalize

IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0);
IRF = IRF_resp / IRF_normalize; % normalize by response of normalization variable
IRF = IRF';

end