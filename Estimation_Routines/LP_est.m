function [IRF,nlags] = LP_est(data_sim,settings);

% preparations

run('Estimation_Setup');

% estimate LP

IRF_resp = IRF_LP(Y,recursiveShock,responseV,nlags,IRF_hor - 1);
IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0);
IRF = IRF_resp / IRF_normalize;
IRF = IRF';

end