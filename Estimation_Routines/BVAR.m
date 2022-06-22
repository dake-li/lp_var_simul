function [Bc,By,Sigma,Beta,Alpha,posteriorVarInv] = BVAR(Y,nlags,prior)
% Auxiliary function for estimating Bayesian VAR coefficients

%% PREPARE

% unpack settings
nv = size(Y,2);
nT = size(Y,1);

towards_random_walk = prior.towards_random_walk;
tight_overall = prior.tight_overall;
tight_nonown_lag = prior.tight_nonown_lag;
decay_power = prior.decay_power;
tight_exogenous = prior.tight_exogenous;

% least-squares VAR to estimate residual variance matrix (treated as known afterwards)
[~,~,Sigma,~,~] = VAR(Y,nlags);
SigmaDiag = diag(Sigma);
SigmaInv = inv(Sigma);

% construct regressors X
X = lagmatrix(Y,1:nlags);
Y = Y((nlags+1):end,:);
X = X((nlags+1):end,:);
X = [ones(size(X,1),1),X];
nx = size(X,2);

%% MN PRIOR (MEAN ZERO / TOWARDS RANDOM WALK)

% prior mean
priorMeanMat = zeros(nv, nv*nlags);
if towards_random_walk == 1
    priorMeanMat(1:nv, 1:nv) = eye(nv);
end

% placeholder for prior variance for each VAR coefficient
priorVarDiagMat = NaN(nv, nv*nlags);

% overall prior tightness
priorVarDiagMat(:) = tight_overall;

% non-own lag prior tightness
nonown_lag_adjustMat = ones(nv, nv) * tight_nonown_lag + eye(nv)*(-tight_nonown_lag+1);
nonown_lag_adjustMat = repmat(nonown_lag_adjustMat, [1, nlags]);
priorVarDiagMat = priorVarDiagMat .* nonown_lag_adjustMat;

% prior variance unit scaling based on residual magnitude
unit_scaling_adjustMat = SigmaDiag * (1./SigmaDiag');
unit_scaling_adjustMat = repmat(unit_scaling_adjustMat, [1, nlags]);
unit_scaling_adjustMat = min(100, max(0.01, unit_scaling_adjustMat)); % truncate by 0.01 and 100 to avoid ill-conditioned matrix
priorVarDiagMat = priorVarDiagMat .* unit_scaling_adjustMat;

% lag decay in prior variance
lag_decay_adjustMat = (1:nlags).^(-decay_power); % exponential decaying with lags
lag_decay_adjustMat = kron(lag_decay_adjustMat, ones(nv));
priorVarDiagMat = priorVarDiagMat .* lag_decay_adjustMat;

% add exogenous variable (constant term) prior mean and tightness
priorMeanMat = [zeros(nv, 1), priorMeanMat];
priorVarDiagMat = [ones(nv,1)*tight_overall*tight_exogenous, priorVarDiagMat];

% construct diagonal of prior variance matrix
priorVarDiag = priorVarDiagMat';
priorVarDiag = priorVarDiag(:);

% construct posterior variance
posteriorVarInv = diag(1./priorVarDiag)+kron(SigmaInv,X'*X);

% compute posterior mean of VAR coefficients alpha = vec(Beta)
Alpha = posteriorVarInv \ (reshape(priorMeanMat', [], 1) ./priorVarDiag + kron(SigmaInv,X)' * Y(:));

% back out VAR coefficients
Beta = reshape(Alpha, [nx, nv]);
Bc = Beta(1,:)'; % constant term
By = reshape(Beta(2:end,:),[nv,nlags,nv]); % lagged term
By = permute(By,[3,1,2]);

end

