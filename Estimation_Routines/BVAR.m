function [Bc,By,Sigma,posteriorVarInv] = BVAR(Y,nlags,prior)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nv = size(Y,2);
nT = size(Y,1);
tight_overall = prior.tight_overall;
tight_nonown_lag = prior.tight_nonown_lag;
decay_power = prior.decay_power;
tight_exogenous = prior.tight_exogenous;

% regular VAR to estimate residual variance
[~,~,Sigma,~,~] = VAR(Y,nlags);
SigmaDiag = diag(Sigma);
SigmaInv = inv(Sigma);

% construct regressors X
X = lagmatrix(Y,1:nlags);
Y = Y((nlags+1):end,:);
X = X((nlags+1):end,:);
X = [ones(size(X,1),1),X];
nx = size(X,2);

% MN prior (mean zero)
priorVarDiagMat = NaN(nv, nv*nlags);

% overall tightness
priorVarDiagMat(:) = tight_overall;

% non-own lag tightness
nonown_lag_adjustMat = ones(nv, nv) * tight_nonown_lag + eye(nv)*(-tight_nonown_lag+1);
nonown_lag_adjustMat = repmat(nonown_lag_adjustMat, [1, nlags]);
priorVarDiagMat = priorVarDiagMat .* nonown_lag_adjustMat;

% unit scaling
unit_scaling_adjustMat = SigmaDiag * (1./SigmaDiag');
unit_scaling_adjustMat = repmat(unit_scaling_adjustMat, [1, nlags]);
unit_scaling_adjustMat = min(100, max(0.01, unit_scaling_adjustMat));
priorVarDiagMat = priorVarDiagMat .* unit_scaling_adjustMat;

% lag decay
lag_decay_adjustMat = (1:nlags).^(-decay_power);
lag_decay_adjustMat = kron(lag_decay_adjustMat, ones(nv));
priorVarDiagMat = priorVarDiagMat .* lag_decay_adjustMat;

% exogenous variable tightness
priorVarDiagMat = [ones(nv,1)*tight_overall*tight_exogenous, priorVarDiagMat];

% construct diagonal of prior variance matrix
priorVarDiag = priorVarDiagMat';
priorVarDiag = priorVarDiag(:);

% construct posterior variance
posteriorVarInv = diag(1./priorVarDiag)+kron(SigmaInv,X'*X);

% compute posterior peak of alpha = vec(Beta)
Alpha = posteriorVarInv \ (kron(SigmaInv,X)' * Y(:));

% back out beta, Bc, By
Beta = reshape(Alpha, [nx, nv]);
Bc = Beta(1,:)';
By = reshape(Beta(2:end,:),[nv,nlags,nv]);
By = permute(By,[3,1,2]);

end

