function [IRF,nlags,Fstat] = SVAR_IV_est(data_sim,settings);

% residualize IV (reg Z on lagged Z and Y)

data_sim = residualize_z_fn(data_sim,settings);

% preparations

run('Estimation_Setup');

% estimate VAR

[VARout,IVout] = SVAR_IV(Y,IV_est_normalize_var_pos,nlags);
ShockVector = IVout.gamma;
IRF = IRF_SVAR(VARout.By,ShockVector,IRF_hor - 1);
IRF = IRF(response_pos,:)';
Fstat = IVout.Fstat_z;

end