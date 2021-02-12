function [Bc,ByCorrect,Sigma,Sxx] = VAR_CorrectBias(Y,nlags)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2); % order (r,x,y,q)
nT = size(Y,1);
[Bc,By,Sigma,Sxx] = VAR(Y,nlags);

VAR_companion_form = zeros(nv * nlags);
VAR_companion_form(1:nv,:) = reshape(-By(:,:,:), [nv, nv*nlags]);
VAR_companion_form((nv+1):end, 1:(nv*(nlags-1))) = eye(nv*(nlags-1));

if abs(eigs(VAR_companion_form,1)) > 1 % if OLS is non-stationary then no correction (Kilian (1998), p. 220)
    
    ByCorrect = By;
    
else

    % standard bias correction
    
    Ay = reshape(By,[nv,nv * nlags]);
    CompAy = kron(diag(ones(1,nlags - 1),-1),eye(nv));
    CompAy(1:nv,:) = Ay;
    E = eig(CompAy);

    sm = zeros(nv * nlags);
    for i = 1:(nv * nlags)
        sm = sm + E(i,1) * inv(eye(nv * nlags) - E(i,1) * CompAy');
    end

    Gamma0 = Sxx(2:end,2:end) - Sxx(2:end,1) * Sxx(2:end,1)';
    G = zeros(nv * nlags);
    G(1:nv,1:nv) = Sigma;

    b = G * (inv(eye(nv * nlags) - CompAy') + CompAy' / (eye(nv * nlags) - CompAy' * CompAy') + sm) / Gamma0;
    Bias = - b / (nT - nlags);

    CompAyCorrect = CompAy - Bias;
    AyCorrect = CompAyCorrect(1:nv,:);
    ByCorrect = reshape(AyCorrect,[nv,nv,nlags]);

    % ensure stationarity of bias-corrected estimate

    indic_stat = 0;

    delta = 1;

    while indic_stat == 0

        VAR_companion_form = zeros(nv * nlags);
        VAR_companion_form(1:nv,:) = reshape(-ByCorrect(:,:,:), [nv, nv*nlags]);
        VAR_companion_form((nv+1):end, 1:(nv*(nlags-1))) = eye(nv*(nlags-1));

        if abs(eigs(VAR_companion_form,1)) < 1
            
            indic_stat = 1;
            
        else
            
            % ad-hoc adjustment as in Kilian (1998), p.220
            
            delta    = delta - 0.01;
            Bias_adj = delta * Bias;

            CompAyCorrect = CompAy - Bias_adj;
            AyCorrect = CompAyCorrect(1:nv,:);
            ByCorrect = reshape(AyCorrect,[nv,nv,nlags]);
            
        end
    end
    
end